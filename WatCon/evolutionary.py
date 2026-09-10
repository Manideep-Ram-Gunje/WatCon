"""Evolutionary conservation (ConSurf) attached to WatCon residues and waters.

This is the join between two things built separately:

* :mod:`WatCon.consurf` parses a ConSurf grades file and can answer
  ``lookup_structural(chain, resid, icode)`` with one of three states.
* :mod:`WatCon.residue_index` gives every WatCon residue an explicit identity,
  ``(chain, resid, icode)``.

Phase 3 joins them on that identity and exposes the result to the network
builders.  It deliberately stops short of inventing a single "conservation
number" for a water or a network -- see *Scientific notes* below.

Scope
-----

This module **consumes ConSurf results that already exist on disk**.  It never
contacts the ConSurf web server.  Submitting a structure and retrieving results
automatically is a separate, later concern; the intended seam for it is
:func:`find_consurf_file`, which is the only place that decides *which* file
belongs to a structure.  A future runtime/API hook should produce a file (or a
parsed result) and hand it to :func:`load_conservation`, leaving everything
below unchanged.

Scientific notes
----------------

**score vs grade.**  ConSurf's ``score`` is already normalised *within a run* to
mean 0 and standard deviation 1 (measured across our four real datasets: mean
-0.000, sd 0.995-0.998).  It is therefore a within-protein z-score of
evolutionary rate -- comparable as a relative rank inside its own protein, and
**not** as an absolute rate between proteins.  ``grade`` is a 1-9 binning of that
score using per-run thresholds, so the bin edges move between runs.

Measured on the same protein run twice at different MSA depths (P00648, 150 vs
50 sequences, 107 shared positions): Spearman 0.955 for score against 0.947 for
grade, and grade changed for 43% of positions.  **Grade is the less stable
representation.**  Both are carried; score is the one to compute with.

**Sign convention.**  More negative score = more conserved.  Grade runs the other
way: 9 = most conserved.  Anything that ranks or colours must not assume "higher
is more conserved".

**Confidence.**  Excluding ConSurf's low-confidence positions barely improved
cross-run agreement in our data (mean grade difference 0.505 -> 0.495), so
filtering costs data and buys little.  Low confidence is recorded on every
residue and counted in aggregates; nothing is dropped here.

**Aggregation.**  A water contacts protein *atoms*, and several atoms may belong
to one residue, so contacts are de-duplicated to residues before aggregating.
No single water-level number is produced -- :class:`WaterConservation` reports
several statistics and lets the consumer choose.  Contact distances are not
available (the network connection tuples do not carry them), so aggregation is
unweighted.

**Naming.**  WatCon already uses "conservation" for *structural water*
conservation (``find_conserved_networks.find_commonality``).  Everything here is
evolutionary conservation and is named ``evo_*`` / ``evolutionary_*`` so the two
never collide in a dict key, a CSV column, or a plot axis.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from statistics import fmean
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

from .consurf import (
    ConSurfParseError,
    ConSurfParseResult,
    ConSurfRecord,
    LookupStatus,
    parse_consurf,
)
from .residue_index import ResidueKey

__all__ = [
    "ConservationError",
    "ResidueConservation",
    "WaterConservation",
    "ConservationMap",
    "ColumnConservation",
    "CROSS_RUN_AGREEMENT",
    "conservation_by_msa_column",
    "family_summary",
    "CoverageReport",
    "DEFAULT_IDENTITY_THRESHOLD",
    "enforce_identity",
    "aggregate_water",
    "aggregate_site",
    "water_residue_contacts",
    "ClusterConservation",
    "conservation_of_clusters",
    "conservation_summary",
    "write_conservation_report",
    "REPORT_COLUMNS",
    "find_consurf_file",
    "load_conservation",
]


#: Fraction of compared residues that must agree before a ConSurf file is
#: accepted for a structure.  Set from measurement, not taste: a correctly
#: numbered structure agrees at 100%, a single point mutant at 99.1% (107/108
#: on barnase), and a one-residue numbering offset collapses to 2.8%.  0.95
#: therefore sits in a wide empty gap -- it tolerates a handful of real
#: mutations while catching any systematic misalignment.
DEFAULT_IDENTITY_THRESHOLD = 0.95


class ConservationError(ValueError):
    """Raised when conservation data cannot be attached safely."""


# ---------------------------------------------------------------------------
# Residue level
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ResidueConservation:
    """ConSurf's verdict for one residue, keyed by structural identity.

    Immutable and shared: every ``OtherAtom`` of a residue points at the *same*
    instance, so the per-atom representation of the network costs one reference
    per atom rather than a copy of six fields.
    """

    #: Normalised conservation score.  MORE NEGATIVE = MORE CONSERVED.
    score: float
    #: ConSurf colour grade, 1-9.  9 = most conserved (opposite sense to score).
    grade: int
    #: ConSurf's own reliability marker (the ``*`` in the COLOR column).
    low_confidence: bool

    #: Non-gapped sequences supporting this alignment column, and the total.
    msa_present: int
    msa_total: int

    #: Bayesian confidence interval on the score (None under maximum likelihood).
    confidence_lower: Optional[float]
    confidence_upper: Optional[float]

    #: ConSurf annotations: 'b'/'e' buried or exposed, 'f'/'s' functional or structural.
    buried_exposed: Optional[str]
    functional_structural: Optional[str]

    #: ConSurf's own POS for this residue.  NOT a dense sequence index --
    #: ConSurf omits non-standard residues, so gaps in POS are expected.
    consurf_position: int

    #: The amino acid ConSurf recorded here, one-letter (its SEQ column).
    #: Used to verify that a structure's residue really is the one this score
    #: describes -- the join key alone cannot detect a numbering offset.
    amino_acid: str = ""
    #: Three-letter name from ConSurf's ATOM column, i.e. what ConSurf actually
    #: read out of *its* PDB.  None for a run with no model attached.
    residue_name: Optional[str] = None
    #: Where the data came from, for provenance.
    source: str = ""

    @property
    def msa_fraction(self) -> float:
        if self.msa_total == 0:
            return 0.0
        return self.msa_present / self.msa_total

    @property
    def is_conserved(self) -> bool:
        """Convenience label for the conserved end (grade 8 or 9).

        Chosen because that is the most reproducible region across ConSurf runs
        in our data (grade >= 8: 26/34 positions identical between two runs of
        the same protein).  This is a label, not a threshold for arithmetic.
        """
        return self.grade >= 8

    @classmethod
    def from_record(cls, record: ConSurfRecord, source: str) -> "ResidueConservation":
        confidence = record.confidence
        return cls(
            score=record.score,
            grade=record.grade,
            low_confidence=record.low_confidence,
            msa_present=record.msa.present,
            msa_total=record.msa.total,
            confidence_lower=None if confidence is None else confidence.score_lower,
            confidence_upper=None if confidence is None else confidence.score_upper,
            buried_exposed=record.buried_exposed,
            functional_structural=record.functional_structural,
            consurf_position=record.position,
            amino_acid=record.sequence_residue,
            residue_name=None if record.pdb_residue is None else record.pdb_residue.name3,
            source=source,
        )

    def to_dict(self) -> dict:
        return {
            "evo_score": self.score,
            "evo_grade": self.grade,
            "evo_low_confidence": self.low_confidence,
            "evo_msa_present": self.msa_present,
            "evo_msa_total": self.msa_total,
            "evo_consurf_position": self.consurf_position,
            "evo_amino_acid": self.amino_acid,
        }


# ---------------------------------------------------------------------------
# Water level
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class WaterConservation:
    """Conservation of the residues one water contacts.

    Several statistics rather than one number, deliberately.  Which of these is
    biologically the right summary is an open question, and baking a single
    choice in now would force a re-run to change it.

    ``n_residues`` counts DISTINCT residues, not contacts: a water touching one
    residue through two atoms scores ``n_residues == 1``.
    """

    #: Most conserved contact (minimum score -- remember, lower = more conserved).
    min_score: float
    #: Mean score over distinct scored contacting residues.
    mean_score: float
    #: Highest grade among contacts (grade 9 = most conserved).
    max_grade: int

    #: Distinct residues with a ConSurf score.
    n_residues: int
    #: How many of those ConSurf flagged low-confidence.
    n_low_confidence: int
    #: Distinct contacting residues with NO ConSurf score (not counted above).
    n_unscored: int

    def to_dict(self) -> dict:
        return {
            "evo_min_score": self.min_score,
            "evo_mean_score": self.mean_score,
            "evo_max_grade": self.max_grade,
            "evo_n_residues": self.n_residues,
            "evo_n_low_confidence": self.n_low_confidence,
            "evo_n_unscored": self.n_unscored,
        }


def aggregate_water(
    conservations: Iterable[Optional[ResidueConservation]],
    n_unscored: Optional[int] = None,
) -> Optional[WaterConservation]:
    """Summarise the residues a water contacts.

    Parameters
    ----------
    conservations:
        One entry per DISTINCT contacting residue.  ``None`` entries are
        residues with no ConSurf score; they are counted, not silently dropped.
        De-duplication is the caller's job -- see
        :meth:`WatCon.generate_static_networks.WaterNetwork.annotate_water_conservation`.
    n_unscored:
        Override for the unscored count.  Normally left as None, in which case
        it is the number of ``None`` entries.

    Returns
    -------
    WaterConservation or None
        None when the water contacts no scored residue at all.  That is a real
        and distinct answer from "contacts conserved residues" -- callers must
        not read it as zero.
    """
    entries = list(conservations)
    scored = [c for c in entries if c is not None]
    unscored = len(entries) - len(scored) if n_unscored is None else n_unscored

    if not scored:
        return None

    scores = [c.score for c in scored]
    return WaterConservation(
        min_score=min(scores),
        mean_score=fmean(scores),
        max_grade=max(c.grade for c in scored),
        n_residues=len(scored),
        n_low_confidence=sum(1 for c in scored if c.low_confidence),
        n_unscored=unscored,
    )


# ---------------------------------------------------------------------------
# The map
# ---------------------------------------------------------------------------

@dataclass
class CoverageReport:
    """How much of a structure a ConSurf run actually covers, and whether the
    residues it covers are the ones the structure actually contains.

    Coverage and identity answer different questions.  Coverage asks *how many*
    residues the run reached; identity asks whether the residue sitting at each
    matched position is the amino acid ConSurf scored there.  A run can cover
    100% of a structure and still be describing a different protein, or the same
    protein numbered differently -- the join key ``(chain, resid, icode)`` cannot
    tell the difference on its own.
    """

    matched: int = 0
    unscored: int = 0
    not_in_file: int = 0
    per_chain: Dict[str, Tuple[int, int]] = None  # chain -> (matched, missing)

    #: Residues (not atoms) whose amino acid was compared and agreed.
    identity_matched: int = 0
    #: Residues whose amino acid disagreed with ConSurf's.
    identity_mismatched: int = 0
    #: (key, structure residue, consurf residue) for each disagreement.
    mismatches: List[Tuple[ResidueKey, str, str]] = None

    def __post_init__(self) -> None:
        if self.per_chain is None:
            self.per_chain = {}
        if self.mismatches is None:
            self.mismatches = []

    @property
    def identity_checked(self) -> int:
        """Residues for which a comparison was actually possible."""
        return self.identity_matched + self.identity_mismatched

    @property
    def identity_rate(self) -> Optional[float]:
        """Fraction of compared residues that agreed, or None if none could be.

        Read this as a rate, not a pass/fail.  A point mutant *should* disagree
        at the mutated position: that is correct data about a real difference
        between the structure and the sequence ConSurf aligned.  A numbering
        offset disagrees almost everywhere.  The two are told apart by where the
        rate falls, which is why this is a number and not a boolean.
        """
        checked = self.identity_checked
        return self.identity_matched / checked if checked else None

    def describe_identity(self) -> str:
        rate = self.identity_rate
        if rate is None:
            return "identity not checked (no residue names available)"
        return (
            f"identity {self.identity_matched}/{self.identity_checked} "
            f"({rate:.1%}) agree"
        )

    @property
    def total(self) -> int:
        return self.matched + self.unscored + self.not_in_file

    @property
    def fraction(self) -> float:
        return self.matched / self.total if self.total else 0.0

    def describe(self) -> str:
        chains = ", ".join(
            f"{c}: {m}/{m + miss}" for c, (m, miss) in sorted(self.per_chain.items())
        )
        return (
            f"matched={self.matched} unscored={self.unscored} "
            f"not_in_file={self.not_in_file} ({self.fraction:.0%} covered)"
            + (f" [{chains}]" if chains else "")
        )


class ConservationMap:
    """Residue identity -> :class:`ResidueConservation`.

    Built from a parsed ConSurf result.  Lookups return ``None`` for a residue
    ConSurf did not score; they never guess.
    """

    def __init__(
        self,
        by_residue: Dict[ResidueKey, ResidueConservation],
        result: Optional[ConSurfParseResult] = None,
        source: str = "",
    ) -> None:
        self._by_residue = dict(by_residue)
        self.result = result
        self.source = source

    # -- construction ------------------------------------------------------

    @classmethod
    def build(
        cls,
        result: ConSurfParseResult,
        chain_map: Optional[Dict[str, str]] = None,
        source: Optional[str] = None,
    ) -> "ConservationMap":
        """Build from a parsed ConSurf result.

        Parameters
        ----------
        result:
            Output of :func:`WatCon.consurf.parse_consurf`.
        chain_map:
            Optional ``{consurf_chain: structure_chain}`` translation, for runs
            whose chain labels differ from the structure being analysed.  The
            default is identity.
        """
        source = source or result.provenance.source
        mapping: Dict[ResidueKey, ResidueConservation] = {}

        for record in result.mapped_records():
            residue = record.pdb_residue
            chain = residue.chain
            if chain_map:
                chain = chain_map.get(chain, chain)
            mapping[(chain, residue.seq, residue.icode)] = (
                ResidueConservation.from_record(record, source)
            )

        return cls(mapping, result=result, source=source)

    # -- lookup ------------------------------------------------------------

    def for_residue(
        self, chain: str, resid: int, icode: Optional[str] = None
    ) -> Optional[ResidueConservation]:
        """Conservation for one residue, or ``None`` when ConSurf did not score it."""
        return self._by_residue.get((chain, resid, icode))

    def for_atom(self, atom) -> Optional[ResidueConservation]:
        """Conservation for the residue an ``OtherAtom`` belongs to."""
        return self.for_residue(
            getattr(atom, "chain", ""),
            atom.resid,
            getattr(atom, "icode", None),
        )

    def status(
        self, chain: str, resid: int, icode: Optional[str] = None
    ) -> LookupStatus:
        """Three-valued answer, distinguishing 'unscored' from 'out of range'.

        ``NO_SCORE_IN_FILE`` means the residue sits inside the range ConSurf
        covered but has no record -- a modified residue ConSurf omitted, such as
        7O7W's chromophore.  ``NOT_IN_FILE`` means it was outside the run
        entirely, such as a chain ConSurf never saw.
        """
        if (chain, resid, icode) in self._by_residue:
            return LookupStatus.SCORED
        if self.result is None:
            return LookupStatus.NOT_IN_FILE
        return self.result.lookup_structural(chain, resid, icode).status

    # -- introspection -----------------------------------------------------

    def __len__(self) -> int:
        return len(self._by_residue)

    def __contains__(self, key: object) -> bool:
        return key in self._by_residue

    def keys(self) -> List[ResidueKey]:
        return list(self._by_residue)

    def chains(self) -> List[str]:
        seen: List[str] = []
        for chain, _, _ in self._by_residue:
            if chain not in seen:
                seen.append(chain)
        return sorted(seen)

    def coverage(self, residues: Sequence) -> CoverageReport:
        """Compare this map against a structure's residues.

        ``residues`` is any sequence of objects exposing ``chain``/``resid``/
        ``icode`` -- :class:`WatCon.residue_index.StructureResidue` fits, and so
        does ``OtherAtom``.  Intended as a pre-flight check: a low coverage
        fraction usually means the ConSurf run and the structure disagree about
        chain labelling.
        """
        report = CoverageReport()
        counts: Dict[str, List[int]] = {}
        # Identity is a property of a residue, but this may be handed one entry
        # per ATOM.  Check each residue once, or a large residue would outvote a
        # small one and the rate would measure atom counts instead of agreement.
        identity_seen: set = set()

        for residue in residues:
            chain = getattr(residue, "chain", "")
            icode = getattr(residue, "icode", None)
            bucket = counts.setdefault(chain, [0, 0])

            status = self.status(chain, residue.resid, icode)
            if status is LookupStatus.SCORED:
                report.matched += 1
                bucket[0] += 1

                key = (chain, int(residue.resid), icode)
                if key not in identity_seen:
                    verdict = self._compare_identity(key, residue)
                    if verdict is not None:
                        identity_seen.add(key)
                        observed, expected = verdict
                        if observed == expected:
                            report.identity_matched += 1
                        else:
                            report.identity_mismatched += 1
                            report.mismatches.append((key, observed, expected))
            else:
                bucket[1] += 1
                if status is LookupStatus.NO_SCORE_IN_FILE:
                    report.unscored += 1
                else:
                    report.not_in_file += 1

        report.per_chain = {c: (m, miss) for c, (m, miss) in counts.items()}
        return report

    def _compare_identity(self, key: ResidueKey, residue) -> Optional[Tuple[str, str]]:
        """(structure residue, ConSurf residue), or None if no comparison is possible.

        Prefers the three-letter name, which is what ConSurf read out of its own
        PDB, and falls back to the one-letter SEQ code for a run with no model
        attached.  Returns None -- rather than guessing or raising -- when the
        object carries no residue name at all, so callers that pass bare
        coordinates degrade to "not checked" instead of "everything mismatched".
        """
        conservation = self._by_residue.get(key)
        if conservation is None:
            return None

        resname = getattr(residue, "resname", None)
        if resname and conservation.residue_name:
            return resname.strip().upper(), conservation.residue_name.strip().upper()

        one_letter = getattr(residue, "one_letter", None)
        if one_letter and conservation.amino_acid:
            return one_letter.strip().upper(), conservation.amino_acid.strip().upper()

        return None


# ---------------------------------------------------------------------------
# Identity enforcement
# ---------------------------------------------------------------------------

def enforce_identity(
    coverage: CoverageReport,
    label: str = "structure",
    strict: bool = True,
    threshold: float = DEFAULT_IDENTITY_THRESHOLD,
) -> Optional[float]:
    """Refuse a ConSurf file whose residues are not the structure's residues.

    The join key ``(chain, resid, icode)`` will happily attach a score to the
    wrong residue if the ConSurf run and the structure disagree about numbering.
    Nothing downstream can detect that: the coverage fraction stays at 100% and
    every water gets a plausible-looking score.  This is the only place the
    mismatch is visible, so it is checked here rather than trusted.

    Returns the identity rate, or None when no comparison was possible (no
    residue names available).  Raises :class:`ConservationError` under ``strict``
    when the rate falls below ``threshold``; warns otherwise.

    A few disagreements are expected and fine -- point mutants genuinely differ
    from the sequence ConSurf aligned.  A systematic offset is not.  See
    :data:`DEFAULT_IDENTITY_THRESHOLD` for why the line sits where it does.
    """
    import warnings

    rate = coverage.identity_rate
    if rate is None or rate >= threshold:
        return rate

    examples = ", ".join(
        "%s%s%s: structure has %s, ConSurf has %s"
        % (key[0], key[1], key[2] or "", observed, expected)
        for key, observed, expected in coverage.mismatches[:5]
    )
    message = (
        "ConSurf data does not describe %s: only %d of %d compared residues "
        "agree (%.1f%%, threshold %.0f%%). This usually means the run and the "
        "structure disagree about residue numbering, in which case every score "
        "attached would be wrong. First disagreements: %s"
        % (label, coverage.identity_matched, coverage.identity_checked,
           rate * 100, threshold * 100, examples or "none recorded")
    )
    if strict:
        raise ConservationError(message)
    warnings.warn(message, RuntimeWarning, stacklevel=2)
    return rate


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

#: Filename suffixes a ConSurf grades file is expected to carry.
CONSURF_SUFFIXES = ("_consurf_grades.txt", ".grades.txt", "_grades.txt")


def find_consurf_file(
    structure_name: str, consurf_directory: Union[str, Path]
) -> Optional[Path]:
    """Locate the ConSurf grades file belonging to one structure.

    Matching uses the same two tokens WatCon's FASTA lookup uses -- the whole
    stem and the part before the first underscore -- but tries the **most
    specific first**, and refuses to choose when a token matches more than one
    file.

    Both of those matter.  Trying the loose token first made ``P00648_50``
    resolve to ``P00648_150.grades.txt``: the bare prefix ``P00648`` matched
    both files and the alphabetically first one won, silently attaching another
    run's conservation to the structure.  Nothing downstream could detect it.

    This function is the single place that decides which ConSurf file belongs to
    a structure, and is therefore the seam where an automatic
    submit-to-ConSurf-and-download step would slot in later: such a hook would
    fetch results into ``consurf_directory`` and leave everything downstream
    untouched.

    Returns
    -------
    Path or None
        None when no candidate matches.  Callers decide whether that is fatal.
    """
    directory = Path(consurf_directory)
    if not directory.is_dir():
        return None

    candidates = [
        p for p in sorted(directory.iterdir())
        if p.is_file() and p.name.endswith(CONSURF_SUFFIXES)
    ]
    if not candidates:
        return None

    stem = Path(structure_name).name
    for suffix in (".pdb", ".cif", ".gro", ".prmtop", ".parm7"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
            break

    # Most specific token first, so a longer name is never beaten by its own
    # prefix matching some other file.
    for token in (stem, stem.split("_")[0]):
        if not token:
            continue
        matches = [
            candidate for candidate in candidates
            if token.lower() in candidate.name.lower()
        ]
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ConservationError(
                "%r matches %d ConSurf files in %s (%s). Rename them so each "
                "structure matches exactly one, rather than having one picked "
                "arbitrarily."
                % (token, len(matches), directory,
                   ", ".join(m.name for m in matches))
            )

    return None


def load_conservation(
    structure_name: str,
    consurf_directory: Optional[Union[str, Path]],
    chain_map: Optional[Dict[str, str]] = None,
    strict: bool = True,
) -> Optional[ConservationMap]:
    """Load the ConSurf results for one structure, if any.

    Parameters
    ----------
    structure_name:
        Structure filename; used to find the matching grades file.
    consurf_directory:
        Directory of ConSurf grades files.  ``None`` disables the whole feature
        and returns ``None``.
    chain_map:
        Optional ``{consurf_chain: structure_chain}`` translation.
    strict:
        ``True`` -- a missing or unparseable file raises
        :class:`ConservationError`.
        ``False`` -- warn and return ``None``, so the run continues with every
        residue unscored.

    Returns
    -------
    ConservationMap or None

    Note
    ----
    Reads supplied files only.  Nothing here contacts the ConSurf server; see
    :func:`find_consurf_file` for where that would attach.
    """
    if consurf_directory is None:
        return None

    path = find_consurf_file(structure_name, consurf_directory)
    if path is None:
        message = (
            f"No ConSurf grades file found for {structure_name!r} in "
            f"{consurf_directory!r}"
        )
        if strict:
            raise ConservationError(message)
        print(f"Warning: {message}; continuing with no conservation data.")
        return None

    try:
        result = parse_consurf(path, strict=strict)
    except ConSurfParseError as exc:
        message = f"Could not parse ConSurf file {path}: {exc}"
        if strict:
            raise ConservationError(message) from exc
        print(f"Warning: {message}; continuing with no conservation data.")
        return None

    return ConservationMap.build(result, chain_map=chain_map, source=str(path))


# ---------------------------------------------------------------------------
# Water -> residue contacts
# ---------------------------------------------------------------------------

def water_residue_contacts(network):
    """Map each water to the DISTINCT residues it contacts.

    Returns ``{water_oxygen_index: {residue_key: ResidueConservation or None}}``.

    De-duplication to residues happens here.  ``OtherAtom`` is a per-atom object,
    so a water touching one residue through both its N and its O would otherwise
    count that residue twice.

    This is the single implementation of the contact walk.  Both network classes'
    ``annotate_water_conservation`` call it, so the two pipelines cannot drift
    apart -- an earlier duplicated version was exactly the kind of thing that
    silently diverges.
    """
    contacts = {}
    connections = getattr(network, "connections", None)
    if not connections:
        return contacts

    atoms_by_index = {atom.index: atom for atom in network.protein_atoms}
    water_by_oxygen = {w.O.index: w for w in network.water_molecules}

    for connection in connections:
        if connection[3] != "WAT-PROT":
            continue

        # Either end may be the protein atom depending on which code path built
        # the connection, so resolve by membership rather than by position.
        first, second = connection[0], connection[1]
        if first in atoms_by_index and second in water_by_oxygen:
            atom, water_index = atoms_by_index[first], second
        elif second in atoms_by_index and first in water_by_oxygen:
            atom, water_index = atoms_by_index[second], first
        else:
            continue

        # int() matters: MDAnalysis hands back numpy integers, which hash like
        # Python ints (so lookups happen to work) but serialise as
        # "np.int64(70)" and would leak numpy types into reports and pickles.
        key = (atom.chain, int(atom.resid), atom.icode)
        contacts.setdefault(water_index, {})[key] = atom.evolutionary

    return contacts


def aggregate_site(network, waters=None):
    """Conservation of the residues lining an arbitrary set of waters.

    Used for the active region, or any hand-picked water selection.  Aggregates
    over DISTINCT residues across all the given waters, so a residue contacted
    by three waters in the site counts once.

    ``waters`` defaults to every water in the network.  Returns ``None`` when no
    contacted residue has a ConSurf score -- "no data", never "not conserved".
    """
    if waters is None:
        waters = network.water_molecules

    contacts = water_residue_contacts(network)
    by_residue = {}
    for water in waters:
        for key, conservation in contacts.get(water.O.index, {}).items():
            # Keep a scored entry in preference to an unscored one for the same
            # residue; they can differ only if the same residue appears twice.
            if key not in by_residue or by_residue[key] is None:
                by_residue[key] = conservation

    return aggregate_water(by_residue.values())


# ---------------------------------------------------------------------------
# Cluster level -- the cross-structure join
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ClusterConservation:
    """Evolutionary conservation of the residues lining one conserved water site.

    A "cluster" here is a WatCon conserved water site: an XYZ position where
    waters recur across a family, produced by
    ``find_conserved_networks.cluster_coordinates_only``.

    This is the object that makes the original scientific question answerable.
    It carries the STRUCTURAL measure (``occupancy``, ``n_structures_occupied``)
    alongside the EVOLUTIONARY one (``min_score`` and friends), side by side and
    never combined -- deciding how they relate is the research question, not
    something the tool should pre-empt.

    Residues are de-duplicated ACROSS structures: if the same residue lines this
    site in five structures, it contributes once.
    """

    cluster_id: int

    # Structural conservation, from WatCon.
    occupancy: int                 # waters at this site, summed over structures
    n_structures_occupied: int     # structures with at least one water here
    n_structures_total: int

    # Evolutionary conservation, from ConSurf.  None when nothing is scored.
    min_score: Optional[float]
    mean_score: Optional[float]
    max_grade: Optional[int]

    n_residues: int                # distinct SCORED residues lining the site
    n_low_confidence: int
    n_unscored: int                # distinct lining residues with no ConSurf score

    residue_keys: Tuple[ResidueKey, ...]

    @property
    def occupancy_fraction(self) -> float:
        """Fraction of structures in which this site is occupied."""
        if self.n_structures_total == 0:
            return 0.0
        return self.n_structures_occupied / self.n_structures_total

    @property
    def has_conservation(self) -> bool:
        return self.min_score is not None

    def to_row(self) -> dict:
        """One flat record, ready for CSV."""
        return {
            "cluster_id": self.cluster_id,
            "occupancy": self.occupancy,
            "n_structures_occupied": self.n_structures_occupied,
            "n_structures_total": self.n_structures_total,
            "occupancy_fraction": round(self.occupancy_fraction, 4),
            "evo_min_score": "NA" if self.min_score is None else round(self.min_score, 4),
            "evo_mean_score": "NA" if self.mean_score is None else round(self.mean_score, 4),
            "evo_max_grade": "NA" if self.max_grade is None else self.max_grade,
            "evo_n_residues": self.n_residues,
            "evo_n_low_confidence": self.n_low_confidence,
            "evo_n_unscored": self.n_unscored,
        }


def _normalise_centers(centers):
    """Accept the several shapes WatCon uses for cluster centres.

    ``cluster_coordinates_only`` returns a ``{label: xyz}`` dict, while
    ``get_coordinates_from_pdb`` returns an ``(N, 3)`` array.  Existing WatCon
    code passes both around interchangeably, so normalise once here to
    ``[(cluster_id, (x, y, z)), ...]`` rather than guessing at each call site.
    """
    if hasattr(centers, "items"):
        return [(int(k), tuple(float(c) for c in v)) for k, v in centers.items()]
    return [(i, tuple(float(c) for c in xyz)) for i, xyz in enumerate(centers)]


def _distance_sq(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def conservation_of_clusters(networks, centers, dist_cutoff=1.5):
    """Join conserved water SITES to the conservation of the residues lining them.

    For every cluster centre, finds the waters within ``dist_cutoff`` in each
    structure, collects the residues those waters contact, de-duplicates them
    across the whole family, and aggregates their ConSurf conservation.

    Parameters
    ----------
    networks : sequence of WaterNetwork
        One per structure.  Conservation must already be attached -- run
        ``initialize_network`` with ``consurf_directory`` set.
    centers : dict or array-like
        Cluster centres from ``find_conserved_networks``.
    dist_cutoff : float
        A water occupies a site if within this distance of the centre, in
        Angstrom.  Matches ``identify_conserved_water_clusters``.

    Returns
    -------
    dict of {cluster_id: ClusterConservation}

    Note
    ----
    A site with no ConSurf-scored lining residue still returns a record, with
    ``min_score=None`` and its ``n_unscored`` count set.  Absence of data is
    reported, never silently rendered as zero conservation.

    Choosing ``dist_cutoff``
    ------------------------
    Cluster centres are the MEAN of their member coordinates, so a loose cluster
    can have a centroid that sits in a gap, further from every water than the
    cutoff.  Observed on real data: clustering one structure's 56 waters gave 15
    centres, 9 of which were more than 1.5 A from any water and so reported
    ``occupancy == 0``.

    That is correct behaviour, not a failure -- but if most sites come back
    unoccupied, the cutoff and the clustering are mismatched rather than the
    protein being interesting.  Cross-family clustering, the intended use,
    produces much tighter sites than clustering a single structure.
    """
    centres = _normalise_centers(centers)
    networks = list(networks)

    # Contact maps are expensive, so build one per structure, not per centre.
    contact_maps = [water_residue_contacts(net) for net in networks]

    # One spatial index per structure, queried once per centre.
    #
    # The obvious triple loop -- every centre against every water of every
    # structure -- is O(centres x structures x waters) in Python. On a real
    # family that is not a small number: 50 PTP1B structures give ~2300 centres
    # and ~11000 waters, which is 27 million distance evaluations for a single
    # join, and the whole point of the tool is to run on families.
    #
    # A KD-tree per structure makes it O(centres x log waters) and returns
    # exactly the same neighbours, so results are unchanged.
    from scipy.spatial import cKDTree

    trees, water_lists = [], []
    for net in networks:
        waters = list(net.water_molecules)
        water_lists.append(waters)
        if waters:
            trees.append(cKDTree([w.O.coordinates for w in waters]))
        else:
            trees.append(None)

    results = {}
    for cluster_id, centre in centres:
        residues = {}
        occupancy = 0
        structures_occupied = 0

        for tree, waters, contacts in zip(trees, water_lists, contact_maps):
            if tree is None:
                continue
            near = tree.query_ball_point(centre, dist_cutoff)
            if not near:
                continue
            occupancy += len(near)
            structures_occupied += 1
            for i in near:
                for key, conservation in contacts.get(waters[i].O.index, {}).items():
                    # Prefer a scored entry: the same residue may be unscored in
                    # one structure and scored in another.
                    if key not in residues or residues[key] is None:
                        residues[key] = conservation

        aggregate = aggregate_water(residues.values())
        unscored = sum(1 for v in residues.values() if v is None)

        results[cluster_id] = ClusterConservation(
            cluster_id=cluster_id,
            occupancy=occupancy,
            n_structures_occupied=structures_occupied,
            n_structures_total=len(networks),
            min_score=None if aggregate is None else aggregate.min_score,
            mean_score=None if aggregate is None else aggregate.mean_score,
            max_grade=None if aggregate is None else aggregate.max_grade,
            n_residues=0 if aggregate is None else aggregate.n_residues,
            n_low_confidence=0 if aggregate is None else aggregate.n_low_confidence,
            n_unscored=unscored,
            residue_keys=tuple(sorted(residues, key=lambda k: (k[0], k[1], k[2] or ""))),
        )

    return results


REPORT_COLUMNS = [
    "cluster_id", "occupancy", "n_structures_occupied", "n_structures_total",
    "occupancy_fraction", "evo_min_score", "evo_mean_score", "evo_max_grade",
    "evo_n_residues", "evo_n_low_confidence", "evo_n_unscored",
]


def conservation_summary(clusters):
    """Cluster records as flat rows, sorted by cluster id."""
    return [clusters[k].to_row() for k in sorted(clusters)]


def write_conservation_report(clusters, path):
    """Write one row per conserved water site.

    Structural conservation (``occupancy``) and evolutionary conservation
    (``evo_*``) sit side by side in the same table, deliberately uncombined.
    Correlating those columns is the scientific question; this file is the input
    to that analysis, not an answer to it.
    """
    import csv

    rows = conservation_summary(clusters)
    directory = os.path.dirname(str(path))
    if directory:
        os.makedirs(directory, exist_ok=True)

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    return len(rows)


# ---------------------------------------------------------------------------
# Family scaffold: pooling conservation from SEPARATE ConSurf runs
# ---------------------------------------------------------------------------
#
# Everything above this line handles one protein: one ConSurf run, applied to
# one or many structures of the same sequence.  That is the barnase case, and it
# needs nothing here.
#
# A true protein FAMILY is different.  Each member has its own ConSurf run, and
# residue 40 of one protein is not residue 40 of another -- the correspondence
# runs through the MSA.  This section provides that join, so adding a family
# means supplying ConSurf files, not writing code.
#
# READ THIS BEFORE POOLING ACROSS RUNS
# ------------------------------------
# ConSurf scores are z-normalised WITHIN each run: across our files the mean is
# -0.000 and the standard deviation 0.995-0.998.  Grades are per-run percentile
# bins of those scores.  Neither is an absolute scale.
#
# So a pooled score answers "how conserved is this position relative to its own
# alignment", never "how conserved is it compared to that other protein".  Two
# runs of the SAME barnase sequence at different MSA depths (150 vs 50
# sequences) agree only at Spearman 0.955 on score and 0.947 on grade, with the
# grade differing at 43% of positions -- and that is the easy case, the same
# sequence.  Different proteins, different alignments, different depths will do
# worse.
#
# The functions below therefore report the spread across members alongside the
# central value, and never average grades.  A grade is an ordinal bin; its mean
# is not a grade.

#: Measured agreement between two ConSurf runs of the same barnase sequence at
#: MSA depths 150 and 50.  Kept here as the honest floor on cross-run pooling:
#: any family-level effect smaller than this is inside ConSurf's own noise.
CROSS_RUN_AGREEMENT = {
    "score_spearman": 0.955,
    "grade_spearman": 0.947,
    "grade_changed_fraction": 0.43,
}


@dataclass(frozen=True)
class ColumnConservation:
    """Conservation at one MSA column, pooled across family members.

    Deliberately not a single number.  ``spread`` is the point: a column where
    every member agrees is evidence, and one where they disagree by more than
    :data:`CROSS_RUN_AGREEMENT` is inside the noise of ConSurf's own
    reproducibility and should not be read as a family signal.
    """

    column: int
    n_members: int
    #: member label -> that member's score at this column
    scores: Dict[str, float]
    #: member label -> that member's grade at this column
    grades: Dict[str, int]
    #: member label -> (chain, resid, icode) contributing the value
    residues: Dict[str, ResidueKey]

    @property
    def mean_score(self) -> Optional[float]:
        return fmean(self.scores.values()) if self.scores else None

    @property
    def min_score(self) -> Optional[float]:
        """Most conserved value across members (score is negative-is-conserved)."""
        return min(self.scores.values()) if self.scores else None

    @property
    def spread(self) -> Optional[float]:
        """Range of scores across members -- how much the runs disagree."""
        if len(self.scores) < 2:
            return None
        return max(self.scores.values()) - min(self.scores.values())

    @property
    def max_grade(self) -> Optional[int]:
        """Highest grade any member assigned.  Grades are never averaged."""
        return max(self.grades.values()) if self.grades else None

    @property
    def unanimous_conserved(self) -> bool:
        """Every member independently called this column conserved (grade >= 8).

        The most defensible family-level statement available, because it does
        not depend on comparing values between separately normalised runs -- only
        on each run's own verdict.
        """
        return bool(self.grades) and all(g >= 8 for g in self.grades.values())


def conservation_by_msa_column(members) -> Dict[int, ColumnConservation]:
    """Pool per-member conservation onto shared MSA columns.

    Parameters
    ----------
    members : sequence of (label, ConservationMap, ResidueIndex)
        One entry per family member.  The ``ResidueIndex`` must carry an MSA
        (``has_msa``), since the alignment column is the only thing that makes
        two different proteins' residues comparable.

    Returns
    -------
    dict of {msa_column: ColumnConservation}

    Raises
    ------
    ConservationError
        If a member has no MSA mapping.  Falling back to residue numbers would
        silently align position 40 of one protein to position 40 of another,
        which is the exact class of error this module exists to prevent.

    Notes
    -----
    Read the section header above before interpreting the output: scores from
    different runs are separately z-normalised and are only relatively
    comparable.  Prefer :attr:`ColumnConservation.unanimous_conserved`, which
    uses each run's own verdict rather than comparing values across runs.
    """
    columns: Dict[int, Dict[str, tuple]] = {}

    for label, conservation_map, index in members:
        if index is None or not index.has_msa:
            raise ConservationError(
                "member %r has no MSA mapping. Residue numbers are not "
                "comparable between different proteins, so pooling without an "
                "alignment would join unrelated positions." % (label,)
            )

        for residue in index.residues:
            column = index.msa_column(residue.chain, residue.resid, residue.icode)
            if column is None:
                continue
            conservation = conservation_map.for_residue(
                residue.chain, residue.resid, residue.icode
            )
            if conservation is None:
                continue
            columns.setdefault(column, {})[label] = (
                conservation.score,
                conservation.grade,
                residue.key,
            )

    result: Dict[int, ColumnConservation] = {}
    for column, per_member in columns.items():
        result[column] = ColumnConservation(
            column=column,
            n_members=len(per_member),
            scores={k: v[0] for k, v in per_member.items()},
            grades={k: v[1] for k, v in per_member.items()},
            residues={k: v[2] for k, v in per_member.items()},
        )
    return result


def family_summary(columns: Dict[int, "ColumnConservation"]) -> dict:
    """Headline counts for a pooled family, including how much runs disagree."""
    if not columns:
        return {
            "n_columns": 0,
            "n_all_members": 0,
            "n_unanimous_conserved": 0,
            "median_spread": None,
        }

    spreads = [c.spread for c in columns.values() if c.spread is not None]
    spreads.sort()
    most = max(c.n_members for c in columns.values())
    return {
        "n_columns": len(columns),
        "n_all_members": sum(1 for c in columns.values() if c.n_members == most),
        "n_unanimous_conserved": sum(
            1 for c in columns.values() if c.unanimous_conserved
        ),
        "median_spread": spreads[len(spreads) // 2] if spreads else None,
    }
