"""Water sites shared across a protein family.

:mod:`WatCon.family` pools *conservation* onto alignment columns. This module
asks the structural question next to it: **does a water recur at the same place
in different proteins of the family, and is that place evolutionarily
conserved?**

Three things have to be right, and each was wrong in an obvious first attempt.

**One frame.** Structures of different proteins are not superposed by residue
number -- ``prepare`` and ``superpose`` both say so, and both assume one protein.
Here the correspondence is the alignment column, so CA atoms of shared columns
are paired and fitted with Kabsch, trimming iteratively to the columns that
actually agree. Measured on the five PTPs against 2F71: the fold core fits at
**0.65-0.90 A** over 178-268 columns, and the active site follows at 0.47-0.93 A
without being fitted directly.

**The WPD loop moves.** In closed structures it sits 0.9-1.3 A from the
reference; in open ones 3.1-8.2 A. That is the conformational change the authors
studied, not a fitting failure -- so a site near that loop means something
different in the two states, and ``per_state_occupancy`` keeps them apart.

**Residues must not merge across proteins.**
:func:`WatCon.evolutionary.conservation_of_clusters` de-duplicates lining
residues by ``(chain, resid, icode)``. Across proteins that is wrong: PTPN1's
Cys215 and PTPN6's residue 215 are unrelated. Here every lining residue is keyed
by ``(protein, alignment column)``, so a site's conservation is stated per
protein and per column, never pooled into one number.

Nothing here invents a combined structural-plus-evolutionary score.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .consurf import parse_consurf
from .evolutionary import ConservationMap, water_residue_contacts
from .family import FamilyConservation, FamilyProtein, FamilyStructure
from .superpose import ca_coordinates, kabsch, rmsd

__all__ = [
    "FAMILY_REPORT_COLUMNS",
    "FamilySiteError",
    "FrameFit",
    "FamilySite",
    "FamilySites",
    "superpose_family",
    "write_superposed",
    "build_family_sites",
    "write_family_report",
]

#: Columns deviating by more than this are dropped before refitting, so a mobile
#: loop cannot drag the whole superposition.
TRIM_CUTOFF = 2.0

#: A fit worse than this over the trimmed core is refused: water sites are
#: clustered at 1.5 A, so a looser frame would invent or destroy sites.
MAX_CORE_RMSD = 1.5

#: Fewest trimmed columns a trustworthy fit may rest on.
MIN_CORE_COLUMNS = 50


class FamilySiteError(ValueError):
    """Raised when structures cannot be brought into one frame."""


@dataclass
class FrameFit:
    """How one structure was placed in the family frame."""

    pdb_id: str
    protein: str
    n_shared_columns: int
    n_core_columns: int
    core_rmsd: float
    rotation: np.ndarray = field(repr=False)
    translation: np.ndarray = field(repr=False)

    def apply(self, coordinates: np.ndarray) -> np.ndarray:
        return np.asarray(coordinates, dtype=float) @ self.rotation.T + self.translation


@dataclass
class FamilySite:
    """One recurring water position, described per protein and per column."""

    site_id: int
    centre: Tuple[float, float, float]

    occupancy: int
    n_structures_occupied: int
    n_structures_total: int
    per_protein_occupancy: Dict[str, int]
    per_state_occupancy: Dict[str, int]
    #: waters at this site contributed by each structure, by PDB id. Kept because
    #: any claim that a site is shared across proteins needs a null, and a null
    #: has to regroup the structures.
    per_structure_occupancy: Dict[str, int]

    #: alignment column -> {protein: best grade of the residue in that column}
    columns: Dict[int, Dict[str, int]]
    #: protein -> the residues lining the site, as (resid, icode)
    residues: Dict[str, Tuple[Tuple[int, Optional[str]], ...]]
    #: lining columns the whole family grades 8 or 9, from the pooled conservation
    unanimous_columns: Tuple[int, ...]
    #: lining residues with no alignment column (excluded conflicts, unplaced)
    n_unplaced_residues: int

    @property
    def proteins_occupied(self) -> Tuple[str, ...]:
        return tuple(sorted(p for p, n in self.per_protein_occupancy.items() if n))

    @property
    def n_proteins_occupied(self) -> int:
        return len(self.proteins_occupied)

    @property
    def max_grade(self) -> Optional[int]:
        grades = [g for by_protein in self.columns.values() for g in by_protein.values()]
        return max(grades) if grades else None

    @property
    def is_family_conserved(self) -> bool:
        """Lined by a column the whole family independently calls conserved."""
        return bool(self.unanimous_columns)


@dataclass
class FamilySites:
    """Every recurring water position in the family, and how it was built."""

    reference: str
    fits: List[FrameFit]
    sites: List[FamilySite]
    n_clusters: int
    n_waters: int
    superposed_dir: Optional[str] = None

    @property
    def shared_sites(self) -> List[FamilySite]:
        """Sites occupied in more than one protein."""
        return [s for s in self.sites if s.n_proteins_occupied > 1]

    def summary(self) -> dict:
        return {
            "n_clusters": self.n_clusters,
            "n_sites_occupied": len(self.sites),
            "n_sites_in_two_or_more_proteins": len(self.shared_sites),
            "n_sites_in_every_protein": sum(
                1 for s in self.sites
                if s.n_proteins_occupied == len({f.protein for f in self.fits})),
            "n_family_conserved_sites": sum(1 for s in self.sites if s.is_family_conserved),
            "n_waters": self.n_waters,
            "worst_core_rmsd": max(f.core_rmsd for f in self.fits) if self.fits else None,
        }


# ---------------------------------------------------------------------------
# One frame
# ---------------------------------------------------------------------------

def _columns_to_ca(family: FamilyConservation, protein: str, structure: FamilyStructure):
    """``{alignment column: CA coordinate}`` for one structure."""
    out: Dict[int, np.ndarray] = {}
    for (_chain, resid, icode), xyz in ca_coordinates(structure.path, chain=structure.chain).items():
        column = family.column_for(protein, resid, icode)
        if column is not None:
            out.setdefault(column, xyz)
    return out


def superpose_family(
    proteins: Sequence[FamilyProtein],
    family: FamilyConservation,
    reference: str,
    trim_cutoff: float = TRIM_CUTOFF,
    max_core_rmsd: float = MAX_CORE_RMSD,
    min_core_columns: int = MIN_CORE_COLUMNS,
) -> List[FrameFit]:
    """Fit every structure onto ``reference`` using CA atoms of shared columns.

    Iteratively drops columns deviating by more than ``trim_cutoff`` and refits,
    so a mobile loop -- the WPD loop moves 3-8 A between open and closed PTPs --
    does not drag the whole structure out of place.

    Raises
    ------
    FamilySiteError
        If the reference is not among the structures, or a structure cannot be
        fitted within ``max_core_rmsd`` over at least ``min_core_columns``.
    """
    by_id = {s.pdb_id: (p.name, s) for p in proteins for s in p.structures}
    if reference not in by_id:
        raise FamilySiteError(
            "reference %r is not one of the structures (%s)"
            % (reference, ", ".join(sorted(by_id))))

    target = _columns_to_ca(family, *by_id[reference])
    if not target:
        raise FamilySiteError("reference %r has no CA atom in any mapped column" % reference)

    fits: List[FrameFit] = []
    for pdb_id, (protein, structure) in sorted(by_id.items()):
        mobile = _columns_to_ca(family, protein, structure)
        shared = sorted(set(mobile) & set(target))
        if len(shared) < min_core_columns:
            raise FamilySiteError(
                "%s shares only %d alignment columns with reference %s (need %d)"
                % (pdb_id, len(shared), reference, min_core_columns))

        keep = shared
        rotation = np.eye(3)
        translation = np.zeros(3)
        for _ in range(5):
            rotation, translation = kabsch(
                np.array([mobile[c] for c in keep]), np.array([target[c] for c in keep]))
            moved = np.array([mobile[c] for c in shared]) @ rotation.T + translation
            deviation = np.linalg.norm(moved - np.array([target[c] for c in shared]), axis=1)
            trimmed = [c for c, d in zip(shared, deviation) if d <= trim_cutoff]
            if len(trimmed) < min_core_columns or trimmed == keep:
                break
            keep = trimmed

        core = rmsd(np.array([mobile[c] for c in keep]) @ rotation.T + translation,
                    np.array([target[c] for c in keep]))
        if core > max_core_rmsd or len(keep) < min_core_columns:
            raise FamilySiteError(
                "%s does not fit the family frame: %.2f A over %d columns "
                "(limits %.2f A, %d columns). Check that it belongs to this family."
                % (pdb_id, core, len(keep), max_core_rmsd, min_core_columns))

        fits.append(FrameFit(pdb_id=pdb_id, protein=protein, n_shared_columns=len(shared),
                             n_core_columns=len(keep), core_rmsd=core,
                             rotation=rotation, translation=translation))
    return fits


def write_superposed(
    proteins: Sequence[FamilyProtein],
    fits: Sequence[FrameFit],
    out_dir: str,
) -> Dict[str, str]:
    """Write each structure transformed into the family frame. Returns paths.

    Coordinates are rewritten in place in the fixed PDB columns; every other
    field, including the chain and any alternate-location label, is left alone.
    """
    by_id = {s.pdb_id: s for p in proteins for s in p.structures}
    os.makedirs(out_dir, exist_ok=True)
    written: Dict[str, str] = {}

    for fit in fits:
        structure = by_id[fit.pdb_id]
        lines = []
        for raw in open(structure.path, "r", errors="replace"):
            line = raw.rstrip("\n")
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                moved = fit.apply(xyz)
                line = (line[:30] + "%8.3f%8.3f%8.3f" % (moved[0], moved[1], moved[2])
                        + line[54:])
            lines.append(line)
        target = os.path.join(out_dir, fit.pdb_id + ".pdb")
        with open(target, "w", newline="\n") as handle:
            handle.write("REMARK   Superposed into the %s family frame: %d columns, RMSD %.2f A.\n"
                         % (fit.protein, fit.n_core_columns, fit.core_rmsd))
            handle.write("\n".join(lines) + "\nEND\n")
        written[fit.pdb_id] = target
    return written


# ---------------------------------------------------------------------------
# The sites
# ---------------------------------------------------------------------------

def build_family_sites(
    proteins: Sequence[FamilyProtein],
    family: FamilyConservation,
    reference: str,
    out_dir: Optional[str] = None,
    states: Optional[Dict[str, str]] = None,
    min_cluster_samples: int = 3,
    site_radius: float = 1.5,
    max_distance: float = 3.3,
    water_name: Optional[str] = "HOH",
    eps: float = 0.0,
    verbose: bool = False,
) -> FamilySites:
    """Superpose the family, cluster its waters, and describe every site.

    Parameters
    ----------
    proteins
        The same members used for :func:`WatCon.family.build_family_conservation`,
        but with ``path`` pointing at structures that still contain their waters.
    family
        The pooled conservation, which supplies the alignment columns.
    reference
        PDB id whose frame everything is placed in.
    states
        Optional ``{pdb id: label}``, e.g. open/closed. Occupancy is reported per
        label, because a site near the WPD loop means different things in the two.
    min_cluster_samples
        Fewest waters forming a site. With one structure per state per protein,
        a site seen in three structures is already cross-protein.
    """
    from scipy.spatial import cKDTree

    from .find_conserved_networks import cluster_coordinates_only
    from .generate_static_networks import extract_objects

    states = dict(states or {})
    fits = superpose_family(proteins, family, reference)
    work = out_dir or os.path.join(os.path.dirname(os.path.abspath(proteins[0].structures[0].path)),
                                   "_family_frame")
    paths = write_superposed(proteins, fits, work)

    protein_of = {f.pdb_id: f.protein for f in fits}
    maps: Dict[str, ConservationMap] = {
        p.name: ConservationMap.build(parse_consurf(p.consurf_path, strict=True)) for p in proteins
    }

    networks, order = [], []
    for fit in fits:
        if verbose:
            print("building network for %s (%s)" % (fit.pdb_id, fit.protein))
        networks.append(extract_objects(
            paths[fit.pdb_id], "water-protein", None,
            active_region_reference=None, active_region_COM=False, active_region_radius=8.0,
            water_name=water_name, msa_indexing=None,
            max_connection_distance=max_distance, conservation_map=maps[fit.protein]))
        order.append(fit.pdb_id)

    coordinates = np.array([w.O.coordinates for net in networks for w in net.water_molecules])
    _labels, centres = cluster_coordinates_only(
        coordinates, cluster="hdbscan", min_samples=min_cluster_samples, eps=eps, source=order)

    contacts = [water_residue_contacts(net) for net in networks]
    waters = [list(net.water_molecules) for net in networks]
    trees = [cKDTree([w.O.coordinates for w in group]) if group else None for group in waters]

    ordered = centres.items() if hasattr(centres, "items") else enumerate(centres)
    sites: List[FamilySite] = []
    for site_id, centre in ordered:
        centre = tuple(float(c) for c in centre)
        occupancy = 0
        occupied_structures = 0
        per_protein: Dict[str, int] = {f.protein: 0 for f in fits}
        per_structure: Dict[str, int] = {}
        per_state: Dict[str, int] = {}
        columns: Dict[int, Dict[str, int]] = {}
        residues: Dict[str, set] = {}
        unplaced = set()

        for index, pdb_id in enumerate(order):
            tree = trees[index]
            if tree is None:
                continue
            near = tree.query_ball_point(centre, site_radius)
            if not near:
                continue
            protein = protein_of[pdb_id]
            occupancy += len(near)
            occupied_structures += 1
            per_protein[protein] += len(near)
            per_structure[pdb_id] = len(near)
            label = states.get(pdb_id)
            if label is not None:
                per_state[label] = per_state.get(label, 0) + len(near)

            for i in near:
                for (_chain, resid, icode), conservation in contacts[index].get(
                        waters[index][i].O.index, {}).items():
                    column = family.column_for(protein, resid, icode)
                    residues.setdefault(protein, set()).add((resid, icode))
                    if column is None:
                        unplaced.add((protein, resid, icode))
                        continue
                    if conservation is None:
                        continue
                    best = columns.setdefault(column, {})
                    if conservation.grade > best.get(protein, 0):
                        best[protein] = conservation.grade

        if occupancy == 0:
            continue

        unanimous = tuple(sorted(
            column for column in columns
            if column in family.columns and family.columns[column].unanimous_conserved))

        sites.append(FamilySite(
            site_id=int(site_id), centre=centre, occupancy=occupancy,
            n_structures_occupied=occupied_structures, n_structures_total=len(order),
            per_protein_occupancy=per_protein, per_state_occupancy=per_state,
            per_structure_occupancy=per_structure,
            columns=columns,
            residues={p: tuple(sorted(v)) for p, v in residues.items()},
            unanimous_columns=unanimous, n_unplaced_residues=len(unplaced),
        ))

    return FamilySites(reference=reference, fits=fits, sites=sites,
                       n_clusters=len(centres), n_waters=len(coordinates),
                       superposed_dir=work)


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

FAMILY_REPORT_COLUMNS = [
    "site_id", "x", "y", "z",
    "occupancy", "n_structures_occupied", "n_structures_total",
    "n_proteins_occupied", "proteins",
    "per_protein_occupancy", "per_structure_occupancy", "per_state_occupancy",
    "n_lining_columns", "n_unanimous_columns", "unanimous_columns",
    "max_grade", "family_conserved", "grades_by_column", "residues_by_protein",
    "n_unplaced_residues",
]


def _pairs(mapping):
    """``a=1;b=2`` -- one CSV cell, sorted, readable in a spreadsheet."""
    return ";".join("%s=%s" % (k, mapping[k]) for k in sorted(mapping))


def write_family_report(sites: FamilySites, path: str) -> int:
    """One row per site. Returns the number of rows written.

    Structural occupancy and evolutionary grades sit in separate columns,
    uncombined, and every grade is attributed to the protein that assigned it --
    the same rule the single-protein report follows.
    """
    import csv

    directory = os.path.dirname(str(path))
    if directory:
        os.makedirs(directory, exist_ok=True)

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FAMILY_REPORT_COLUMNS)
        writer.writeheader()
        for site in sorted(sites.sites, key=lambda s: (-s.n_proteins_occupied, -s.occupancy)):
            grades = {column: _pairs(by_protein) for column, by_protein in site.columns.items()}
            writer.writerow({
                "site_id": site.site_id,
                "x": round(site.centre[0], 3),
                "y": round(site.centre[1], 3),
                "z": round(site.centre[2], 3),
                "occupancy": site.occupancy,
                "n_structures_occupied": site.n_structures_occupied,
                "n_structures_total": site.n_structures_total,
                "n_proteins_occupied": site.n_proteins_occupied,
                "proteins": ";".join(site.proteins_occupied),
                "per_protein_occupancy": _pairs(site.per_protein_occupancy),
                "per_structure_occupancy": _pairs(site.per_structure_occupancy),
                "per_state_occupancy": _pairs(site.per_state_occupancy),
                "n_lining_columns": len(site.columns),
                "n_unanimous_columns": len(site.unanimous_columns),
                "unanimous_columns": ";".join(str(c) for c in site.unanimous_columns),
                "max_grade": "NA" if site.max_grade is None else site.max_grade,
                "family_conserved": "yes" if site.is_family_conserved else "no",
                "grades_by_column": _pairs(grades),
                "residues_by_protein": _pairs({
                    protein: "+".join(str(resid) for resid, _icode in residues)
                    for protein, residues in site.residues.items()}),
                "n_unplaced_residues": site.n_unplaced_residues,
            })
    return len(sites.sites)
