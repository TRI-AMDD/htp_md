"""
Implemented functions to compute ion clusters' population matrix either from Pizza package when aggregates.txt is
available or from trajectory dump file otherwise.
"""
import os
import sys
import os.path
import numpy as np
import itertools
from htpmd.shared.pizza_dump import dump
sys.path.append(".")
sys.path.append("../")


def compute_pop_matrix_from_pizza_dump(trajectory_reference_path, li_idx=90, o_idx=94):
    """
    Description: Receives the path to simulation folder and computes population matrix using dump pizza_dump module

        This function assumes Li ("90") and N in TFSI ("93") if the meta.json file does not exist

    Args: trajectory_reference_path: path to the trajectory file to be processed (Ex:
    '0623/poly_9_nion_50_conc_1.5/9-0-485080326-0')

    Returns:
        pop_matrix: population matrix that includes the number of clusters with different shapes.
        pop_matrix[i,j] is the number of clusters with i number cations and j number of anions
    """

    path_to_aggregate_file = os.path.join(trajectory_reference_path, 'aggregates.txt')

    # Cutoff in the cluster population: don't output a cluster matrix larger than (cutoff,cutoff)
    cutoff = 50

    print('Loading file...')
    f = dump(path_to_aggregate_file)
    print('Done')

    n_snap = f.nsnaps
    tsteps = f.time()

    # Determine max cluster size based on the number of ions
    ids = f.vecs(tsteps[0], 'id')

    f.tselect.none()

    n_max = int(np.size(ids) / 2)
    popmatrix = np.zeros((max(cutoff, n_max), max(cutoff, n_max)), dtype=float)

    # Loop over all snapshots, create intermediate list with
    # unique cluster ids, then cast list into the population matrix
    for i in range(0, n_snap):
        f.tselect.one(tsteps[i])
        time = f.time()[0]

        ccids = np.unique(f.vecs(time, 'c_cc2'))
        atoms = np.asarray(f.vecs(time, 'type', 'c_cc2', 'mol'))

        for j in range(0, np.size(ccids)):
            catccids = atoms[1, np.where(atoms[0, :] == li_idx)]
            aniccids = atoms[1, np.where(atoms[0, :] == o_idx)]

            catmols = atoms[2, np.where(atoms[0, :] == li_idx)]
            animols = atoms[2, np.where(atoms[0, :] == o_idx)]

            ccmols = catmols[np.where(catccids == ccids[j])]
            acmols = animols[np.where(aniccids == ccids[j])]

            ncat = len(np.unique(ccmols))
            nani = len(np.unique(acmols))

            popmatrix[ncat, nani] += 1.0

    # We now have the whole population matrix. Output as i j pop
    popmatrix /= n_snap

    # Apply cutoff
    popmatrix = popmatrix[:cutoff, :cutoff]

    return popmatrix


def add_info(info1, info2):
    """
    Description:
    Receives two pieces of information (can be tuples or lists,), concatenate, and return the joined info as a tuple
    This is used to add values calculated for the current timestep to the previous results.

    Args:
        info1 (tuple or list)
        info2 (tuple or list)

    Returns:
        joined_info a tuple which is concatenated form of info1 and info2

    """
    joined_info = None

    if isinstance(info1, list) and isinstance(info2, list):
        joined_info = info1 + info2

    elif isinstance(info1, tuple) and isinstance(info2, list):
        joined_info = list(info1) + info2

    elif isinstance(info1, list) and isinstance(info2, tuple):
        joined_info = info1 + list(info2)

    elif isinstance(info1, tuple) and isinstance(info2, tuple):
        joined_info = list(info1) + list(info2)

    else:
        print('Unexpected data!')
    return tuple(joined_info)


def distance(vec1, vec2):
    """
    Description:
    Receives coordinates of two points and computes the distance between two of them (ends of two vectors)

    Args:
        vec1 (coordinates of point 1 with respect to the origin)
        vec2 (coordinates of point 2 with respect to the origin)

    Returns:
        d (distance between point 1 and point 2)

    """

    d = np.sqrt((vec1[0] - vec2[0]) ** 2
                + (vec1[1] - vec2[1]) ** 2
                + (vec1[2] - vec2[2]) ** 2)

    return d


def wrap_atom_coords(step_data, box, atom_id):
    """
    Description:
        Receives step_data that includes coordinates of atoms, the box info including the size of the box, and atom
        id and wrap the coordinate.
        We noticed that in some cases the wrapped coordinates obtained from LAMMPS does not work as expected
        and the atoms are still outside of the box.

    Args:
        step_data - tuples of information extracted from the current time step
        box - includes the size of the box in the x, y, and z directions
        atom_id
    Returns:
        x, y, z (wrapped coordinates)
    """

    box_x = box[0]['max'] - box[0]['min']
    box_y = box[1]['max'] - box[1]['min']
    box_z = box[2]['max'] - box[2]['min']

    # wrap atom coordinates into the simulation box
    if step_data[atom_id][5] > box[0]['max']:
        x = step_data[atom_id][5] - box_x
    elif step_data[atom_id][5] < box[0]['min']:
        x = step_data[atom_id][5] + box_x
    else:
        x = step_data[atom_id][5]

    if step_data[atom_id][6] > box[1]['max']:
        y = step_data[atom_id][6] - box_y
    elif step_data[atom_id][6] < box[1]['min']:
        y = step_data[atom_id][6] + box_y
    else:
        y = step_data[atom_id][6]

    if step_data[atom_id][7] > box[2]['max']:
        z = step_data[atom_id][7] - box_z
    elif step_data[atom_id][7] < box[2]['min']:
        z = step_data[atom_id][7] + box_z
    else:
        z = step_data[atom_id][7]
        #
    return x, y, z


def shortest_distance(step_data, box, atom_id1, atom_id2):
    """
    Description:
        Received step_data that includes coordinates of atoms, box info that includes the size of the box, and atom id
        and atom_id2 and returns the shortest distance between atoms with id1 and id2 in the periodic box

    Args:
        step_data: tuples of information extracted from the current time step
        box: box parameters including min and max position of x, y, and z
        atom_id1: atom of interest for which we are looking for neighbors
        atom_id2: can be the neighbor of atom_id1 or not, based on the shortest distance


    Returns:
        the shortest distance between atom_id1 and atom_id2 in a periodic box
    """
    periodic_distances = []

    box_x = box[0]['max'] - box[0]['min']
    box_y = box[1]['max'] - box[1]['min']
    box_z = box[2]['max'] - box[2]['min']

    x_1, y_1, z_1 = wrap_atom_coords(step_data, box, atom_id1)
    x_2, y_2, z_2 = wrap_atom_coords(step_data, box, atom_id2)

    wall_id_array = np.zeros((1, 3))

    # identify the wall proximity in x, y, and z directions

    if x_2 >= box_x / 2:
        wall_id_array[0, 0] = -1
    else:
        wall_id_array[0, 0] = 1

    if y_2 >= box_y / 2:
        wall_id_array[0, 1] = -1
    else:
        wall_id_array[0, 1] = 1

    if z_2 >= box_z / 2:
        wall_id_array[0, 2] = -1
    else:
        wall_id_array[0, 2] = 1

    x_combination = [x_2, x_2 + wall_id_array[0, 0] * box_x]

    y_combination = [y_2, y_2 + wall_id_array[0, 1] * box_y]

    z_combination = [z_2, z_2 + wall_id_array[0, 2] * box_z]

    for p_image in list(itertools.product(*[x_combination, y_combination, z_combination])):
        d = distance([x_1, y_1, z_1], list(p_image))
        if d <= box_x / 2:
            return d
        else:
            periodic_distances.append(d)

    # Sort distances from smallest to largest
    periodic_distances.sort()

    return periodic_distances[0]


def get_cluster_info(step_data, box, trace_ids, trace_ids2):
    """
    Description: Receives atom data for the current time step(a tuple), and two lists of atom ids (anions and
    cations), and identifies the clusters formed from combination of two lists - trace_ids("Li 90") and trace_ids2 (
    atom IDs of O, S, and N in TFSI)  based on a cut-off radius of 3.4 Angstrom.

        for more info, please see the following reference:

        France-Lanord and Grossman. "Correlations from ion pairing and the
        Nernst-Einstein equation." Physical review letters 122.13 (2019): 136001.

    Args:
        step_data: tuple of atoms' information extracted from the current time step

        box: box parameters including min and max position of x, y, and z

        trace_ids: List of atom Ids with the type "90" (Li atoms)

        trace_ids2: List of atom Ids with the type "93", "94", or "95" (O, S, and N in TFSI)
        we needed to have all these to that avoid counting TFSI ions several times. (distance between O atoms at two
        ends > 3.4 Angstrom)

    Returns:
        all_clusters: a list of lists. individual lists are clusters identified in the current timestep
    """

    merged_list = trace_ids + trace_ids2
    all_clusters = []

    while merged_list != []:
        this_cluster = []
        this_cluster.append(merged_list[0])
        merged_list.remove(merged_list[0])

        for i in this_cluster:
            for j in merged_list:
                d = shortest_distance(step_data, box, i, j)
                if d <= 3.4:
                    this_cluster.append(j)
                    merged_list.remove(j)

        all_clusters.append(this_cluster)

    return all_clusters


def population_matrix(step_data, all_clusters, type_id=90, type_id3=93):
    """
    Description: Receives atom data for the current time step(a tuple), atom ids that formed clusters in this time
    step( all_clusters: a list of lists)

        and type of atoms for Li ("90") and N in TFSI ("93")

        I used N in TFSI because there is only one N in each TFSI

    Args:
        step_data: tuple of atoms' information extracted from the current time step

        all_clusters: a list of lists. individual lists are clusters identified in the current timestep

        type_id: "90" (Li atoms)

        type_id3: "93" (N in TFSI)

    Returns:
        pop_matrix: population matrix that includes the number of clusters with different shapes.
        pop_matrix[i,j] is the number of clusters with i number cations and j number of anions
    """

    # Li : type_id = 90 , N: type_id3 = 93
    pop_matrix = np.zeros((50, 50, 1))

    for cluster in all_clusters:
        cations = 0
        anions = 0

        for atom_id in cluster:

            if step_data[atom_id][2] == type_id:
                cations += 1
            if step_data[atom_id][2] == type_id3:
                anions += 1

        pop_matrix[cations][anions] += 1

    return pop_matrix


def stack_population_matrix(stacked_population, current_population):
    """
    Description: Takes a stack of clusters' population matrices (a NumPy array) and updates it by adding the current
    population matrice to that.

    Args:
        stacked_population: a stack of population matrices from previous time frames
        current_population: current population matrix related to the current time frame
    Returns:
        stacked_population: updated population matrices including the data from the current time frame
    """
    if not np.size(stacked_population):
        stacked_population = current_population

    else:
        stacked_population = np.dstack((stacked_population, current_population))

    return stacked_population


def generate_population_matrix(trajectory_reference_path, type_id=90, type_id2=94, type_id3=95, type_id4=93,
                               min_steps=0, max_steps=2500, step_size=1000):
    """
    Description:

        Iterates over the timesteps in the trajectory and computes mean squared displacements and population
         matrices.

        Args:
        trajectory_reference: path to the trajectory file to be processed
             Ex: '0623/poly_9_nion_50_conc_1.5/9-0-485080326-0'

        type_id: "90" (Li atoms)

        type_id2: "94" (O in TFSI)

        type_id3: "93" (N in TFSI)

        type_id4: "95" (S in TFSI)

        n_nearest: the number of nearest neighbors to track (currently not used)

        n_samples: the number of sampled Li ions to track (currently not used)

        min_steps: starting time frame number to gather information

        max_steps: final time frame number to gather information

        step_size: Frequency of output in trajectory or multiplication of frequency by an integer

    Returns:
        Saves the averaged population matrix from the beginning up to the end
        and the stack of population matrices for each time step
    """
    trajectory_file_path = f'{trajectory_reference_path}/traj.lammpstrj'
    trajectory_trace_path = f'{trajectory_reference_path}/traj_trace.data'

    dirname = os.path.dirname(trajectory_trace_path)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    f_target = open(trajectory_trace_path, 'w')

    f_target.write('id,mol,type,mass,q,x,y,z,ix,iy,iz,x_unfold,y_unfold,z_unfold,timestep\n')

    # initialize 3D population matrix for all time steps
    stacked_population = np.array([])
    #

    with open(trajectory_file_path, 'r') as f:
        data_lines = f.readlines()
        tot_lines = len(data_lines)

        # steps through the data file and separates the data into time step groups
        # assuming same sized header block for each time step
        l_num = 0  # line number
        while l_num < tot_lines:

            line = data_lines[l_num]
            if line[:14] == 'ITEM: TIMESTEP':  # starting a new time step

                trace_ids = []  # initializing the list of trace particles - Li
                trace_ids2 = []  # initializing the list of trace particles - O in TFSI

                step_num = int(data_lines[l_num + 1])

                if (step_num >= min_steps * step_size) and (step_num <= max_steps * step_size) and (
                        step_num % step_size == 0):

                    # gets the box size
                    box = {
                        0: {
                            'min': None,
                            'max': None
                        },
                        1: {
                            'min': None,
                            'max': None
                        },
                        2: {
                            'min': None,
                            'max': None
                        }
                    }

                    size_x = (data_lines[l_num + 5].split())
                    size_y = (data_lines[l_num + 6].split())
                    size_z = (data_lines[l_num + 7].split())

                    box[0]['min'], box[0]['max'] = [float(x) for x in size_x]
                    box[1]['min'], box[1]['max'] = [float(x) for x in size_y]
                    box[2]['min'], box[2]['max'] = [float(x) for x in size_z]
                    ###

                    step_data = {}

                    num_atoms = int(data_lines[l_num + 3])

                    l_num = l_num + 9

                    # grab atom data for this time step
                    for atom in range(num_atoms):
                        line = data_lines[l_num + atom].split()
                        line_data = tuple(map(float, line[:]))
                        atom_id = line_data[0]
                        step_data[atom_id] = line_data[0:]
                        line_type_id = line_data[2]

                        # identify Li ions
                        if line_type_id == type_id:
                            trace_ids.append(atom_id)

                        # identify O, N and S in TFSI
                        if line_type_id == type_id2 or line_type_id == type_id3 or line_type_id == type_id4:
                            trace_ids2.append(atom_id)

                        # compute and add the unwarpped coordinates to the step_data
                        x = line_data[5]
                        ix = line_data[8]
                        x_unfold = box[0]['min'] + ix * (box[0]['max'] - box[0]['min']) + x
                        y = line_data[6]
                        iy = line_data[9]
                        y_unfold = box[1]['min'] + iy * (box[1]['max'] - box[1]['min']) + y
                        z = line_data[7]
                        iz = line_data[10]
                        z_unfold = box[2]['min'] + iz * (box[2]['max'] - box[2]['min']) + z

                        step_data[atom_id] = add_info(step_data[atom_id], [x_unfold, y_unfold, z_unfold])

                    # Makes a tuple from information of Li ions and O, S, and N in TFSI and step number
                    for this_id in trace_ids + trace_ids2:
                        step_data[this_id] = add_info(step_data[this_id], [step_num])
                        f_target.write(','.join(map(str, step_data[this_id])) + '\n')

                    # Trace clusters
                    all_clusters = get_cluster_info(step_data, box, trace_ids, trace_ids2)

                    # form population matrix using clusters' info
                    current_population = population_matrix(step_data, all_clusters, type_id=90, type_id3=93)

                    # form stacked population matrix using population matrix at each time frame
                    stacked_population = stack_population_matrix(stacked_population, current_population)

                else:

                    l_num = l_num + 1

            else:
                l_num = l_num + 1

        # save population matrix averaged from min_steps to max_steps
        avg_population = np.mean(stacked_population, axis=2)

        return stacked_population, avg_population
