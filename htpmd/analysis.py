import json
import os
import numpy as np
from .utils import load_lammps, _ATOM_MASSES


REQUIRED_METADATA = {
    'mol_smiles',  # Smiles for molecule. E.g. '[Cu]CCO[Au]' for PEO
    'poly_smiles',  # Smiles for polymer.  E.g. 'CCOCCOCCOCCOCCOCCO' for PEO
    'force_field',  # Force field. E.g. 'PCFF+'
    'material_group',  # Type of materials. E.g. polymer
    'temperature',  # Temperature in K. E.g. 353
    'timestep',  # timestep in fs. E.g. 2.
}


def get_metadata(dir_name):
    with open(os.path.join(dir_name, 'meta.json')) as f:
        metadata = json.load(f)
    assert set(metadata.keys()).issuperset(REQUIRED_METADATA)
    return metadata


def get_diffusivity(dir_name, target_type=90):
    """Diffusivity of a specified atom type (unit: cm^2/s)."""
    _, _, types, unwrapped_coords = load_lammps(
        os.path.join(dir_name, 'traj.lammpstrj'))
    return compute_diff(types, unwrapped_coords, target_type)


def get_molarity(dir_name):
    """Molarity of the polymer/salt mixture (unit: mol Li / kg polymer)."""
    lammps_file = os.path.join(dir_name, 'traj.lammpstrj')
    _, _, raw_types, _ = load_lammps(lammps_file, use_mass=False)
    _, _, real_types, _ = load_lammps(lammps_file, use_mass=True)

    poly_idx = np.nonzero(raw_types < 89)[0]
    li_idx = np.nonzero(real_types == 3)[0]

    atom_masses = np.array(_ATOM_MASSES)

    poly_mass = np.sum(atom_masses[real_types[poly_idx]])

    return float(len(li_idx)) / poly_mass * 1e3


def compute_diff(types, unwrapped_coords, target_type):
    target_idx = np.nonzero(types == target_type)[0]
    target_coords = unwrapped_coords[:, target_idx]
    msd = np.mean(np.sum((target_coords[-1] - target_coords[0])**2, axis=-1))
    return msd / (len(target_coords) - 1) / 6 * 5e-5 # cm^2/s


def get_conductivity(dir_name):
    """Total conducitivty of the system (unit: S/cm)."""
    max_cluster = 10
    e_const = 1.6021766209e-19
    kb_const = 1.38064852e-23
    T = 353.0

    pop_mat = np.loadtxt(os.path.join(dir_name, 'population.txt'))

    _, lattices, types, unwrapped_coords = load_lammps(
        os.path.join(dir_name, 'traj.lammpstrj'))

    li_diff = compute_diff(types, unwrapped_coords, target_type=90)  # cm^2/s
    tfsi_diff = compute_diff(types, unwrapped_coords, target_type=93)  # cm^2/s

    assert np.isclose(lattices[0:1], lattices).all()

    V = np.prod(lattices[0]) * 1e-24  # cm^3

    cond = 0.
    total_ion = 0.

    for i in range(max_cluster):
        for j in range(max_cluster):
            if i > j:
                cond += e_const**2 / V / kb_const / T * \
                    (i - j)**2 * pop_mat[i, j] * li_diff
            elif i < j:
                cond += e_const**2 / V / kb_const / T * \
                    (i - j)**2 * pop_mat[i, j] * tfsi_diff
            else:
                pass
            total_ion += (i + j) * pop_mat[i, j]

    return cond  # S/cm


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
