"""Dump each residue's secondary-structure assignment from PyMOL.

A small standalone helper, inherited from upstream and not imported by
the rest of the package. Run it inside PyMOL, with a structure loaded::

    run print_structure_from_pymol.py

It writes one line per residue -- name, number, and PyMOL's ``ss``
assignment (H helix, S sheet, L loop) -- which is useful for reading
alongside a conservation or occupancy table.
"""

from pymol import cmd

def save_secondary_structure(outfile="ss_output.txt", selection="name CA"):
    """Write the secondary structure of a selection to a text file.

    Parameters
    ----------
    outfile : str, optional
        Where to write. One line per residue.
    selection : str, optional
        PyMOL selection to iterate. The default picks one atom per
        residue, which is what makes the output one line per residue
        rather than one per atom.
    """
    with open(outfile, "w") as f:
        def writer(resn, resi, ss):
            f.write(f"{resn}{resi:>4}  {ss}\n")
        cmd.iterate(selection, "writer(resn, resi, ss)", space={"writer": writer})

if __name__ == '__main__':
    save_secondary_structure()