"""Water sites shared across the PTP family, from real structures.

Five proteins, one closed structure each, trimmed to 12 A around their own
catalytic nucleophile with the waters kept, plus the CA extracts and the
published alignment that place them in one frame.

The claim under test is narrow and checkable: after superposing on CA atoms of
shared alignment columns, waters from five different proteins land in the same
places, and the places lined by residues every ConSurf run calls conserved
include the catalytic water -- the position lined by the Asp181, Cys215 and
Gln262 columns.

The subtle failure these tests guard against is residues merging across
proteins. ``conservation_of_clusters`` de-duplicates lining residues by
``(chain, resid, icode)``; across proteins that would join PTP1B's Asp181 to
whatever PTPN6 numbers 181. Here every site keeps each protein's own residue
numbers, and only alignment columns are shared.
"""

from __future__ import annotations

import os

import pytest

pytest.importorskip("Bio", reason="family mapping needs Biopython")
pytest.importorskip("MDAnalysis", reason="building networks needs MDAnalysis")

from WatCon.family import FamilyProtein, FamilyStructure, build_family_conservation
from WatCon.family_sites import (
    FamilySiteError,
    build_family_sites,
    superpose_family,
    write_superposed,
)

from .conftest import ptp_grades

PACKAGE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAMILY_DIR = os.path.join(PACKAGE, "data", "examples", "ptp_family")
ALIGNMENT = os.path.join(FAMILY_DIR, "ptp_family_alignment.pir")

#: protein -> (ConSurf run, closed structure, its trimmed active-site fixture)
SITES = {
    "PTPN1": ("1AAX", "2F71", "2F71_A_site.pdb"),
    "PTPN6": ("4GRZ", "4GRZ", "4GRZ_A_site.pdb"),
    "PTPN7": ("1ZC0", "1ZC0", "1ZC0_A_site.pdb"),
    "PTPN12": ("5HDE", "5HDE", "5HDE_A_csp_site.pdb"),
    "PTPN22": ("3BRH", "3BRH", "3BRH_A_site.pdb"),
}
#: the CA extracts used only to build the alignment columns
PAIRS = {
    "PTPN1": ("2F71", "8U1E"), "PTPN6": ("4GRZ", "4HJP"), "PTPN7": ("1ZC0", "3O4U"),
    "PTPN12": ("5HDE", "5J8R"), "PTPN22": ("3BRH", "3OLR"),
}

pytestmark = pytest.mark.skipif(
    not all(os.path.isfile(os.path.join(FAMILY_DIR, f)) for _r, _c, f in SITES.values()),
    reason="PTP family site fixtures absent")


@pytest.fixture(scope="module")
def family():
    """Conservation and columns, from the full-chain CA extracts."""
    return build_family_conservation(
        [FamilyProtein(name, str(ptp_grades(SITES[name][0])),
                       [FamilyStructure(pdb_id, os.path.join(FAMILY_DIR, "%s_A_ca.pdb" % pdb_id))
                        for pdb_id in pair],
                       reference=SITES[name][1])
         for name, pair in PAIRS.items()],
        ALIGNMENT)


@pytest.fixture(scope="module")
def site_proteins():
    """The trimmed active-site structures, which carry the waters."""
    return [FamilyProtein(name, str(ptp_grades(run)),
                          [FamilyStructure(closed, os.path.join(FAMILY_DIR, fixture))],
                          reference=closed)
            for name, (run, closed, fixture) in SITES.items()]


@pytest.fixture(scope="module")
def sites(family, site_proteins, tmp_path_factory):
    out = tmp_path_factory.mktemp("family_frame")
    return build_family_sites(site_proteins, family, reference="2F71",
                              out_dir=str(out), min_cluster_samples=2)


# ===========================================================================
# One frame
# ===========================================================================

def test_five_proteins_fit_one_frame(family, site_proteins):
    fits = superpose_family(site_proteins, family, reference="2F71", min_core_columns=30)
    assert len(fits) == 5
    assert all(fit.n_core_columns >= 55 for fit in fits)
    assert all(fit.core_rmsd < 0.65 for fit in fits)


def test_the_reference_is_not_moved(family, site_proteins):
    fits = superpose_family(site_proteins, family, reference="2F71", min_core_columns=30)
    reference = next(f for f in fits if f.pdb_id == "2F71")
    assert reference.core_rmsd == pytest.approx(0.0, abs=1e-9)


def test_a_fit_that_cannot_be_trusted_is_refused(family, site_proteins):
    """Demanding more columns than the trimmed structures share must stop, not guess."""
    with pytest.raises(FamilySiteError, match="alignment columns"):
        superpose_family(site_proteins, family, reference="2F71", min_core_columns=1000)


def test_superposed_files_keep_every_atom_and_move_them(family, site_proteins, tmp_path):
    fits = superpose_family(site_proteins, family, reference="2F71", min_core_columns=30)
    written = write_superposed(site_proteins, fits, str(tmp_path))

    def atoms(path):
        return [l for l in open(path) if l.startswith(("ATOM", "HETATM"))]

    for fit in fits:
        source = next(s for p in site_proteins for s in p.structures if s.pdb_id == fit.pdb_id)
        before, after = atoms(source.path), atoms(written[fit.pdb_id])
        assert len(before) == len(after), fit.pdb_id
        moved = sum(1 for a, b in zip(before, after) if a[30:54] != b[30:54])
        if fit.pdb_id == "2F71":
            assert moved == 0                      # the reference stays put
        else:
            assert moved == len(before), fit.pdb_id


# ===========================================================================
# The sites
# ===========================================================================

def test_the_family_has_sites_every_protein_occupies(sites):
    summary = sites.summary()
    assert summary["n_waters"] == 81
    assert summary["n_sites_in_every_protein"] >= 5
    assert summary["worst_core_rmsd"] < 0.65


def _catalytic_columns(family):
    return {family.column_for("PTPN1", resid) for resid in (181, 215, 262)}


def test_the_catalytic_water_recurs_in_all_five_proteins(sites, family):
    """The site lined by the Asp181, Cys215 and Gln262 columns, in every protein."""
    wanted = _catalytic_columns(family)
    catalytic = [s for s in sites.sites
                 if wanted <= set(s.columns) and s.n_proteins_occupied == 5]
    assert catalytic, "no site lined by all three catalytic columns in all five proteins"
    site = catalytic[0]
    assert set(site.proteins_occupied) == set(SITES)
    assert site.unanimous_columns
    assert site.is_family_conserved


#: One alignment column, five different residue numbers -- each protein's own
#: equivalent. Corroborated independently: the nucleophile column lands on the
#: residue each PDB entry declares (PTPN6 C453S, PTPN22 C227S), and the Asp181
#: column lands on PTPN22's engineered D195A position.
COLUMN_EQUIVALENTS = {
    181: {"PTPN1": 181, "PTPN6": 419, "PTPN7": 236, "PTPN12": 199, "PTPN22": 195},
    215: {"PTPN1": 215, "PTPN6": 453, "PTPN7": 270, "PTPN12": 231, "PTPN22": 227},
    262: {"PTPN1": 262, "PTPN6": 500, "PTPN7": 314, "PTPN12": 278, "PTPN22": 274},
}


@pytest.mark.parametrize("ptp1b_resid", sorted(COLUMN_EQUIVALENTS))
def test_a_column_means_a_different_residue_in_each_protein(family, ptp1b_resid):
    """The merge this design exists to prevent: 181 in PTP1B is not 181 elsewhere."""
    column = family.column_for("PTPN1", ptp1b_resid)
    observed = {protein: key[1] for protein, key in family.columns[column].residues.items()}
    assert observed == COLUMN_EQUIVALENTS[ptp1b_resid]


def test_each_protein_keeps_its_own_residue_numbers(sites, family):
    """A site's lining residues stay in each protein's numbering, never pooled."""
    wanted = _catalytic_columns(family)
    site = next(s for s in sites.sites if wanted <= set(s.columns) and s.n_proteins_occupied == 5)

    numbers = {protein: {resid for resid, _icode in residues}
               for protein, residues in site.residues.items()}
    assert len(numbers) == 5
    # PTPN6 is numbered in its full-length protein's frame, the others are not.
    assert min(numbers["PTPN6"]) > 400
    assert all(max(numbers[p]) < 400 for p in ("PTPN1", "PTPN7", "PTPN12", "PTPN22"))
    # No two proteins contribute the same set, which a merge would produce.
    as_tuples = [tuple(sorted(v)) for v in numbers.values()]
    assert len(set(as_tuples)) == len(as_tuples)


def test_occupancy_adds_up(sites):
    for site in sites.sites:
        assert sum(site.per_structure_occupancy.values()) == site.occupancy
        assert sum(site.per_protein_occupancy.values()) == site.occupancy
        assert site.n_structures_occupied == len(site.per_structure_occupancy)


def test_states_are_reported_separately(family, site_proteins, tmp_path):
    """Open and closed mean different things near the WPD loop, so they are split."""
    states = {pdb_id: "closed" for _r, pdb_id, _f in SITES.values()}
    built = build_family_sites(site_proteins, family, reference="2F71",
                               out_dir=str(tmp_path), min_cluster_samples=2, states=states)
    assert all(set(s.per_state_occupancy) <= {"closed"} for s in built.sites)
    assert sum(s.per_state_occupancy.get("closed", 0) for s in built.sites) > 0


def test_conserved_sites_are_occupied_by_more_proteins(sites):
    """The family-level finding, in miniature on the fixtures.

    On the full ten-structure dataset: sites lined by a unanimously conserved
    column are occupied by 3.43 proteins on average against 2.66 for the rest
    (Mann-Whitney p = 7.6e-8; 20% vs 3% occupied in all five). Burial is not
    disentangled -- conserved residues are also more buried.
    """
    scored = [s for s in sites.sites if s.columns]
    conserved = [s for s in scored if s.unanimous_columns]
    other = [s for s in scored if not s.unanimous_columns]
    assert conserved and other
    mean = lambda group: sum(s.n_proteins_occupied for s in group) / len(group)
    assert mean(conserved) > mean(other)
