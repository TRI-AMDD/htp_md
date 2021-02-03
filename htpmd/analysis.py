import json
import os
import numpy as np
from .utils import load_lammps
from pymatgen.core.structure import Structure
from htpmd.constants import ATOM_MASSES
from htpmd.trajectory.load import LammpsTrajectoryLoader
from htpmd.shared.polymer import (
    compute_diffusivity, compute_polymer_diffusivity, compute_molarity,
    compute_conductivity, compute_msd_curve, get_cif_at_frame)


REQUIRED_METADATA = {
    'mol_smiles',  # Smiles for molecule. E.g. '[Cu]CCO[Au]' for PEO
    'poly_smiles',  # Smiles for polymer.  E.g. 'CCOCCOCCOCCOCCOCCO' for PEO
    'force_field',  # Force field. E.g. 'PCFF+'
    'material_group',  # Type of materials. E.g. polymer
    'temperature',  # Temperature in K. E.g. 353
    'time_step',  # timestep in fs. E.g. 2.
}


def get_all_properties(dir_name):
    # Load trajectories and pre-computed properties
    traj = LammpsTrajectoryLoader().load(dir_name)
    pop_mat = get_population_matrix(dir_name)
    metadata = get_metadata(dir_name)

    # Compute properties
    results = dict()
    results.update(metadata)
    traj.remove_drift()
    results['li_diffusivity'] = compute_diffusivity(traj, target_type=90)
    results['tfsi_diffusivity'] = compute_diffusivity(traj, target_type=93)
    results['poly_diffusivity'] = compute_polymer_diffusivity(traj)
    results['conductivity'] = compute_conductivity(traj, pop_mat=pop_mat)
    results['molarity'] = compute_molarity(traj)
    results['li_msd_curve'] = compute_msd_curve(traj, target_type=90)
    results['tfsi_msd_curve'] = compute_msd_curve(traj, target_type=93)
    results['structure'] = get_cif_at_frame(traj, k=0)
    return results


def get_metadata(dir_name):
    with open(os.path.join(dir_name, 'meta.json')) as f:
        metadata = json.load(f)
    assert set(metadata.keys()).issuperset(REQUIRED_METADATA)
    return metadata


def get_population_matrix(dir_name):
    """Load the population matrix computed by lammps."""
    pop_mat = np.loadtxt(os.path.join(dir_name, 'population.txt'))
    return pop_mat


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
