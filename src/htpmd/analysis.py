import os
import numpy as np
from .utils import load_lammps
from htpmd.constants import RawType
from htpmd.trajectory.load import (
    LammpsTrajectoryLoader, get_metadata, get_population_matrix)
from htpmd.shared.polymer import (
    compute_diffusivity, compute_polymer_diffusivity, compute_molality,
    compute_conductivity, compute_msd_curve, compute_non_avg_msd_curve, get_cif_at_frame,
    compute_displacement, compute_simulation_length)
from htpmd.ml_models import gnn, random_forest


ML_PROPERTIES = [
    'conductivity',
    'li_diffusivity',
    'poly_diffusivity',
    'tfsi_diffusivity',
    'molarity',
    'transference_number',
]


def get_all_properties(dir_name):
    # Load trajectories and pre-computed properties
    traj = LammpsTrajectoryLoader().load(dir_name)
    pop_mat = get_population_matrix(dir_name)
    metadata = get_metadata(dir_name)

    # Get parameters from metadata
    cation_raw_type = metadata['cation_raw_type']
    anion_raw_type = metadata['anion_raw_type']

    # Compute properties
    results = dict()
    results.update(metadata)
    traj.remove_drift()
    results['li_diffusivity'] = compute_diffusivity(traj, target_type=cation_raw_type, **metadata)
    results['tfsi_diffusivity'] = compute_diffusivity(traj, target_type=anion_raw_type, **metadata)
    results['poly_diffusivity'] = compute_polymer_diffusivity(traj, **metadata)
    results['conductivity'], results['transference_number'] = compute_conductivity(
        traj, pop_mat=pop_mat, **metadata)
    results['molality'] = compute_molality(traj, **metadata)
    results['li_msd_curve'] = compute_msd_curve(traj, target_type=cation_raw_type, **metadata)
    results['tfsi_msd_curve'] = compute_msd_curve(traj, target_type=anion_raw_type, **metadata)
    results['li_msd_curve_non_avg'] = compute_non_avg_msd_curve(traj, target_type=cation_raw_type, **metadata)
    results['tfsi_msd_curve_non_avg'] = compute_non_avg_msd_curve(traj, target_type=anion_raw_type, **metadata)
    results['structure'] = get_cif_at_frame(traj, k=0)
    results['li_mean_disp'] = compute_displacement(traj, target_type=cation_raw_type, type='mean')
    results['li_max_disp'] = compute_displacement(traj, target_type=cation_raw_type, type='max')
    results['tfsi_mean_disp'] = compute_displacement(traj, target_type=anion_raw_type, type='mean')
    results['tfsi_max_disp'] = compute_displacement(traj, target_type=anion_raw_type, type='max')
    results['simulation_length'] = compute_simulation_length(traj, **metadata)

    # Only predict properties for polymers
    if metadata['material_group'] == 'polymer':
        if metadata['mol_smiles'] is not None:
            # Get GNN predicted properties
            for prop in ML_PROPERTIES:
                gnn_pred = gnn.predict([metadata['mol_smiles']], prop)[0]
                results[f'gnn_{prop}'] = gnn_pred
            # Get RF predicted properties
            for prop in ML_PROPERTIES:
                if prop == 'molarity':
                    continue
                rf_pred = random_forest.random_forests_prediction([metadata['mol_smiles']], prop)[0]
                results[f'rf_{prop}'] = rf_pred
    
    return results


def lammpstraj2npz(dir_name, out_dir, target_atom_num):
    """Converts a lammpstraj file to a npz file."""
    lammps_file = os.path.join(dir_name, 'traj.lammpstrj')
    output_file = os.path.join(out_dir, dir_name.split('/')[-1])

    print('Loading lammps file.')
    print('-' * 80)
    traj_coords, lattices, atom_types, unwrapped_coords = load_lammps(
        lammps_file, use_mass=True, tol=0.01)
    target_index = np.nonzero(atom_types == target_atom_num)[0]

    data_dict = {
        'traj_coords': traj_coords,
        'lattices': lattices,
        'atom_types': atom_types,
        'target_index': target_index,
        'unwrapped_coords': unwrapped_coords,
    }
    print('Saving data.')
    print('-' * 80)
    np.savez_compressed(output_file, **data_dict)


def analyze_all(root_dir, analyze_fn, num_workers=None):
    """Function to analyze all directories in"""
    dir_names = [fn for fn in os.listdir(root_dir)
                 if os.path.isdir(os.path.join(root_dir, fn))]

    results = {}

    if num_workers is None:
        for fn in dir_names:
            dir_name = os.path.join(root_dir, fn)
            result = analyze_fn(dir_name)
            results[fn] = result
            print(fn, 'done')
    else:
        from multiprocessing import Process, Manager

        def f(results, fns):
            for fn in fns:
                dir_name = os.path.join(root_dir, fn)
                try:
                    result = analyze_fn(dir_name)
                    results[fn] = result
                    print(fn, 'done')
                except OSError:
                    print(fn, 'OSError')

        with Manager() as manager:
            results_m = manager.dict()

            splited_fnames = [fns.tolist()
                              for fns in np.array_split(dir_names, num_workers)]
            procs = []

            for fns in splited_fnames:
                procs.append(Process(target=f, args=(results_m, fns)))

            for p in procs:
                p.start()
            for p in procs:
                p.join()

            results = dict(results_m)

    return results
