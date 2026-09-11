'''
Create pdbs and pml files for visualization of water networks
'''

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def pymol_project_oxygen_network(network, filename='STATE.pml', out_path='pymol_projections', active_region_only=False, water_only=False):
    """
    Generate a PyMOL (.pml) file to visualize the calculated network.

    Parameters
    ----------
    network : WaterNetwork
        The WaterNetwork object representing the network.
    filename : str
        Name of the output .pml file.
    out_path : str
        Directory where the output file will be saved.
    active_region_only : bool, optional
        Whether to include only active site connections. Default is False.
    water_only : bool, optional
        Whether to show only water-water connections. Default is False.

    Returns
    -------
    None
    """

    os.makedirs(out_path, exist_ok=True)
    with open(os.path.join(out_path, filename), 'w') as FILE:
        #for i, mol in enumerate(network.water_molecules):
            #FILE.write(f"show spheres, id {mol.O.index+1}\nset sphere_scale, 0.3, id {mol.O.index+1}\n")

        #If active_region_only, then only project connections within the active site. Note this is redundant if you have 
        #chosen to only include active site atoms in your whole network
        if active_region_only:
            connection_list = [f for f in network.connections if f[4]=='active_region']
        else:
            connection_list = network.connections

        if water_only:
            connection_list = [f for f in connection_list if f[3]=='WAT-WAT']


        for i, connection in enumerate(connection_list):
            FILE.write(f"distance interaction{i}, id {connection[0]}, id {connection[1]}\n")
            FILE.write(f"show spheres, id {connection[0]}\nshow spheres, id {connection[1]}\n")

        #Stylistic preferences
        FILE.write(f"set dash_radius, 0.15, interaction*\nset dash_color, black, interaction*\n")
        FILE.write(f"set dash_gap, 0.0, interaction*\nhide labels, interaction*\n")
        FILE.write(f"set sphere_scale, 0.2\n")
        FILE.write(f"bg white\n")

        #Create WaterNetwork group
        FILE.write(f"group WaterNetwork, interaction*\n")


def project_clusters(coordinate_list, filename_base='CLUSTER',b_factors=None):
    """
    Generate a PDB file to visualize cluster centers.

    Parameters
    ----------
    coordinate_list : array-like, dict
        List of coordinates representing cluster centers.
    filename_base : str
        Naming scheme to use for the outputted PDB file.
    b_factors : array-like, optional
        Optional list of values to replace the B-factor column.

    Returns
    -------
    None
    """
    
    #Obtain coordinates from cluster centers
    if type(coordinate_list) == dict:
        coordinate_list = coordinate_list.values() 

    if b_factors is None:
        b_factors = [0.00] * len(coordinate_list)  # Default B-factor is 0.00 if not provided

    os.makedirs('cluster_pdbs', exist_ok=True)

    filename = f"cluster_pdbs/{filename_base}.pdb"
    with open(filename, 'w') as FILE:
        atom_serial = 1
        for label, (center, b_factor) in enumerate(zip(coordinate_list, b_factors)):
            # Write each cluster center as an oxygen atom in PDB format
            FILE.write(
                f"ATOM  {atom_serial:5d}  O   HOH A{label:4d}    "
                # PDB fixed columns: occupancy 55-60, tempFactor 61-66.
                f"{center[0]:8.3f}{center[1]:8.3f}{center[2]:8.3f}{1.00:6.2f}{b_factor:6.2f}           O\n"
            )
            atom_serial += 1


def plot_consevation_angles(cluster_conservation_dict, output_filebase='angle_clusters', output_dir='pymol_projections'):
    """
    Visualize two-angle water-protein conservation

    Parameters
    ----------
    cluster_conservation_dict : dict
        Dictionary obtained from WatCon.find_conserved_networks.identify_clustered_angles()
    output_filebase : str, optional
        Filename base (no extension) for outputted .pml file. Default is 'angle_clusters'
    output_dir : str, optional
        Name of output directory. Default is 'pymol_projections'

    Returns
    -------
    None
    """
    os.makedirs(output_dir, exist_ok=True)

    strengths = []


    with open(os.path.join(output_dir, f"{output_filebase}.pml"), 'w') as FILE:
        for msa in cluster_conservation_dict.keys():
            for i, label in enumerate(cluster_conservation_dict[msa].keys()):
                strength = cluster_conservation_dict[msa][label]['counts']
                strengths.append(int(strength))

    strengths = np.array(strengths)

    cmap = plt.get_cmap('bwr')
    len_colors = np.linspace(0,1,len(strengths))
    colors = cmap(len_colors)

    #colors = [(int(r*255), int(g*255), int(b*255)) for r, g, b, _ in colors]
    colors = [(float(r), float(g), float(b)) for r, g, b, _ in colors]
    colors = [color for _, color in sorted(zip(strengths, colors), key=lambda x: x[0])]

    tick = 0
    with open(os.path.join(output_dir, f"{output_filebase}.pml"), 'w') as FILE:
        for msa in cluster_conservation_dict.keys():
            for i, label in enumerate(cluster_conservation_dict[msa].keys()):
                water_coord = cluster_conservation_dict[msa][label]['wat_coord']
                prot_coord = cluster_conservation_dict[msa][label]['closest_coord']

                strength = cluster_conservation_dict[msa][label]['counts']

                FILE.write(f"pseudoatom {msa}_{i}_protein, pos=[{prot_coord[0]}, {prot_coord[1]}, {prot_coord[2]}]\n")
                FILE.write(f"pseudoatom {msa}_{i}_water, pos=[{water_coord[0]}, {water_coord[1]}, {water_coord[2]}]\n")
                FILE.write(f"distance interaction_{msa}_{i}, {msa}_{i}_protein, {msa}_{i}_water\n")
                FILE.write(f"set dash_color, [{colors[tick][0]},{colors[tick][1]},{colors[tick][2]}], interaction_{msa}_{i}\n")
                tick += 1

        
        FILE.write("hide labels, all\n")
        FILE.write('set dash_radius, 0.15, interaction*\n')    
        FILE.write('set dash_gap, 0.0, interaction*\n')
        FILE.write('hide labels, interaction*\n')
        FILE.write('set sphere_scale, 0.4\n')
        FILE.write('set sphere_color, oxygen, *_water')
        FILE.write('bg white\n')
        FILE.write('group AngleInteractions, interaction*\n')
        FILE.write('group PseudoProteins, *_protein\n')
        FILE.write('group PseudoWater, *_water\n')
        FILE.write('show spheres, PseudoProteins or PseudoWater\n')
        FILE.write('color oxygen, PseudoWater\n')

def export_graph_to_pdb(graph, output_file):
    """
    Generate a PDB file with dummy oxygen atoms at the coordinates of a given graph.

    Parameters
    ----------
    graph : networkx.Graph
        The NetworkX graph containing node coordinates.
    output_file : str
        Name of the output PDB file.

    Returns
    -------
    None
    """

    if not output_file.endswith('.pdb'):
        output_file = f"{output_file}.pdb"

    os.makedirs('graph_pdbs', exist_ok=True)
    with open(f"graph_pdbs/{output_file}", 'w') as f:
        atom_serial = 1
        # Collect nodes involved in active site edges
        active_region_nodes = set()
        for edge1, edge2, data in graph.edges(data=True):
            if data.get('active_region') == 'active_region':
                active_region_nodes.add(edge1)
                active_region_nodes.add(edge2)

        # Write nodes as PDB atoms
        for node, data in graph.nodes(data=True):
            if node in active_region_nodes:  # Check if node is part of active site
                x, y, z = data.get('pos', (0.0, 0.0, 0.0))  # Default position if not provided
                f.write(
                    f"ATOM  {atom_serial:5d}  O   HOH A{atom_serial:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O\n"
                )
                atom_serial += 1

def save_xyz(network, filename='STATE.xyz'):
    """
    Generate an XYZ file containing graph coordinates.

    Parameters
    ----------
    network : WaterNetwork
        The WaterNetwork object containing graph data.
    filename : str
        Name of the output XYZ file.

    Returns
    -------
    None
    """
    if not output_file.endswith('.xyz'):
        output_file = f"{output_file}.xyz"

    with open(filename, 'w') as FILE:
        FILE.write(f"{len(network.molecules)}\n")
        for i, mol in enumerate(network.molecules):
            x, y, z  = mol.coordinates
            FILE.write(f"{mol.name}\t{x:.3f} {y:.3f} {z:.3f}\n")



# ---------------------------------------------------------------------------
# Evolutionary conservation projections
# ---------------------------------------------------------------------------
#
# WatCon already projects STRUCTURAL water conservation into the B-factor column
# (see project_clusters).  These functions do the same for EVOLUTIONARY
# conservation, into SEPARATE files.  The two are never blended into one number:
# whether they agree is the scientific question, so a user should be able to
# load both and look.

#: ConSurf's own colour scale, grade 1 (variable, cyan) to 9 (conserved, maroon).
#: Reproduced so projections match what users see on the ConSurf server.
CONSURF_GRADE_COLORS = {
    1: (0.06, 0.78, 0.82),
    2: (0.56, 0.94, 0.94),
    3: (0.85, 1.00, 1.00),
    4: (0.94, 0.94, 1.00),
    5: (1.00, 1.00, 1.00),
    6: (0.98, 0.86, 0.90),
    7: (0.98, 0.63, 0.75),
    8: (0.94, 0.35, 0.58),
    9: (0.63, 0.16, 0.37),
}


def project_clusters_by_conservation(clusters, centers, filename_base='EVO_CLUSTER',
                                     out_dir='cluster_pdbs', value='grade'):
    """Write cluster centres as a PDB with EVOLUTIONARY conservation in B-factors.

    Mirrors :func:`project_clusters`, which carries *structural* conservation in
    the same column.  Writing a separate file keeps the two measures distinct --
    load both in PyMOL and compare.

    Parameters
    ----------
    clusters : dict of {cluster_id: ClusterConservation}
        From ``WatCon.evolutionary.conservation_of_clusters``.
    centers : dict or array-like
        The same centres used to build ``clusters``.
    filename_base : str
        Output basename; ".pdb" is appended.
    out_dir : str
        Output directory.
    value : {'grade', 'score'}
        Which measure goes into the B-factor column.  'grade' is 1-9 with 9 most
        conserved and reads naturally in a viewer.  'score' is the continuous
        measure, where MORE NEGATIVE means MORE conserved -- colour ramps must be
        inverted for it.

    Returns
    -------
    str
        Path written.

    Note
    ----
    Sites with no conservation data get B-factor 0.00 and are listed in the
    REMARK block.  Zero here means "no data", NOT "not conserved"; the remark
    exists so that distinction survives into the file.
    """
    import os

    if value not in ('grade', 'score'):
        raise ValueError("value must be 'grade' or 'score'")

    if hasattr(centers, 'items'):
        ordered = [(int(k), v) for k, v in centers.items()]
    else:
        ordered = list(enumerate(centers))

    os.makedirs(out_dir, exist_ok=True)
    filename = os.path.join(out_dir, filename_base + '.pdb')

    unscored = []
    with open(filename, 'w') as FILE:
        FILE.write('REMARK   B-factor column carries EVOLUTIONARY conservation (%s).\n'
                   % value)
        if value == 'score':
            FILE.write('REMARK   MORE NEGATIVE = MORE CONSERVED.\n')
        else:
            FILE.write('REMARK   Grade 9 = most conserved, 1 = most variable.\n')
        FILE.write('REMARK   B-factor 0.00 means NO DATA, not "not conserved".\n')

        serial = 1
        for cluster_id, center in ordered:
            record = clusters.get(cluster_id)
            if record is None or not record.has_conservation:
                b_factor = 0.00
                unscored.append(cluster_id)
            elif value == 'grade':
                b_factor = float(record.max_grade)
            else:
                b_factor = float(record.min_score)

            FILE.write(
                'ATOM  %5d  O   HOH A%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           O\n'
                % (serial, cluster_id, center[0], center[1], center[2], 1.00, b_factor)
            )
            serial += 1

        if unscored:
            FILE.write('REMARK   sites with no conservation data: %s\n'
                       % ', '.join(str(c) for c in unscored))

    return filename


def pymol_project_evolutionary(network, filename='EVO_RESIDUES.pml',
                               out_path='pymol_projections', grade_cutoff=None,
                               structure=None):
    """Write a .pml colouring protein residues by ConSurf grade.

    Uses ConSurf's own 1-9 colour scale so the projection matches what the user
    sees on the ConSurf server.

    Parameters
    ----------
    network : WaterNetwork
        Built with a conservation map attached.
    filename : str
        Output .pml name.
    out_path : str
        Output directory.
    grade_cutoff : int, optional
        If given, only colour residues at or above this grade -- useful for
        highlighting just the conserved end, which is also the most reproducible
        part of a ConSurf run.
    structure : str, optional
        Path to the structure these residues belong to.  **Strongly recommended.**
        Without it the script is a bare list of ``color`` commands that assume a
        structure is already open; opened on its own it colours an empty session
        and shows nothing at all:

            $ pymol EVO_RESIDUES.pml
            objects loaded: []
            atoms visible: 0

        Every command succeeds, so there is no error to notice.  Passing the
        structure makes the file self-contained.  ``watcon view`` does this for
        you.

    Returns
    -------
    str
        Path written.  Residues with no ConSurf score are left uncoloured and
        reported in a comment, never coloured as if they were variable.
    """
    import os

    os.makedirs(out_path, exist_ok=True)
    filename = os.path.join(out_path, filename)

    # One entry per residue, not per atom: OtherAtom is a per-atom object.
    by_residue = {}
    unscored = set()
    for atom in network.protein_atoms:
        key = (atom.chain, int(atom.resid))
        if atom.evolutionary is None:
            unscored.add(key)
        else:
            by_residue[key] = atom.evolutionary.grade

    by_grade = {}
    for (chain, resid), grade in by_residue.items():
        if grade_cutoff is not None and grade < grade_cutoff:
            continue
        by_grade.setdefault(grade, []).append((chain, resid))

    with open(filename, 'w') as FILE:
        FILE.write('# Evolutionary conservation (ConSurf grade), 9 = most conserved.\n')
        FILE.write('# Residues with no ConSurf score are NOT coloured.\n')
        if structure is not None:
            FILE.write('load %s\n' % os.path.abspath(structure).replace('\\', '/'))
        else:
            FILE.write('#\n')
            FILE.write('# NOTE: no structure was given, so this file only COLOURS.\n')
            FILE.write('# Load your structure first, or it will colour an empty\n')
            FILE.write('# session and show nothing:   pymol structure.pdb this.pml\n')
        FILE.write('bg white\n')
        FILE.write('hide everything\n')
        FILE.write('show cartoon\n')
        FILE.write('color grey80\n')

        for grade in sorted(by_grade):
            r, g, b = CONSURF_GRADE_COLORS[grade]
            FILE.write('set_color consurf_%d, [%.3f, %.3f, %.3f]\n' % (grade, r, g, b))
            for chain, resid in sorted(by_grade[grade]):
                chain_sel = ' and chain %s' % chain if chain.strip() else ''
                FILE.write('color consurf_%d, resi %d%s\n' % (grade, resid, chain_sel))

        FILE.write('# residues coloured: %d\n' % sum(len(v) for v in by_grade.values()))
        FILE.write('# residues with no ConSurf score: %d\n' % len(unscored))

    return filename
