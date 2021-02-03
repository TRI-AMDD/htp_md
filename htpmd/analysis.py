import json
import os
import numpy as np
from .utils import load_lammps
from pymatgen.core.structure import Structure
from htpmd.constants import ATOM_MASSES
from htpmd.trajectory.load import LammpsTrajectoryLoader
from htpmd.shared.transport import compute_diffusivity, compute_polymer_diffusivity


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
    results['conductivity'] = get_conductivity(
        traj.lattices, traj.raw_types, traj.unwrapped_coords, pop_mat)
    results['molarity'] = get_molarity(traj.raw_types, traj.atom_types)
    results['li_msd_curve'] = get_msd_curve(
        traj.raw_types, traj.unwrapped_coords, target_type=90)
    results['tfsi_msd_curve'] = get_msd_curve(
        traj.raw_types, traj.unwrapped_coords, target_type=93)
    results['structure'] = get_cif_at_frame(
        traj.wrapped_coords, traj.lattices, traj.atom_types, k=0)
    return results


def get_metadata(dir_name):
    with open(os.path.join(dir_name, 'meta.json')) as f:
        metadata = json.load(f)
    assert set(metadata.keys()).issuperset(REQUIRED_METADATA)
    return metadata


def get_raw_traj(dir_name):
    """Load the entire trajectory into memory."""
    wrapped_coords, lattices, types, atom_types, unwrapped_coords = load_lammps(
        os.path.join(dir_name, 'traj.lammpstrj'))
    return wrapped_coords, lattices, types, atom_types, unwrapped_coords


def get_population_matrix(dir_name):
    """Load the population matrix computed by lammps."""
    pop_mat = np.loadtxt(os.path.join(dir_name, 'population.txt'))
    return pop_mat


def get_molarity(raw_types, real_types):
    """Molarity of the polymer/salt mixture (unit: mol Li / kg polymer)."""

    poly_idx = np.nonzero(raw_types < 89)[0]
    li_idx = np.nonzero(real_types == 3)[0]

    atom_masses = np.array(ATOM_MASSES)

    poly_mass = np.sum(atom_masses[real_types[poly_idx]])

    return float(len(li_idx)) / poly_mass * 1e3


def compute_center_of_mass(coords, atom_types):
    element_masses = np.array(ATOM_MASSES)
    atom_masses = element_masses[atom_types]
    return (np.sum(coords * atom_masses[np.newaxis, :, np.newaxis], axis=1) /
            np.sum(atom_masses))


def get_coords_without_drift(coords, atom_types):
    cos_coord = compute_center_of_mass(coords, atom_types)
    return coords - cos_coord[:, np.newaxis]


def get_diffusivity(raw_types, unwrapped_coords, target_type):
    """Diffusivity of a specified atom type (unit: cm^2/s)."""
    target_idx = np.nonzero(raw_types == target_type)[0]
    target_coords = unwrapped_coords[:, target_idx]
    msd = np.mean(np.sum((target_coords[-1] - target_coords[0])**2, axis=-1))
    return msd / (len(target_coords) - 1) / 6 * 5e-5  # cm^2/s


def get_polymer_diffusivity(raw_types, atom_types, unwrapped_coords):
    # TODO: should F be included?
    solvate_types = (atom_types == 7) | (atom_types == 8) | (atom_types == 16)
    poly_solvate_types = (raw_types < 90) & solvate_types
    poly_solvate_idx = np.nonzero(poly_solvate_types)[0]
    target_coords = unwrapped_coords[:, poly_solvate_idx]
    msd = np.mean(np.sum((target_coords[-1] - target_coords[0])**2, axis=-1))
    return msd / (len(target_coords) - 1) / 6 * 5e-5  # cm^2/s


def get_msd_curve(raw_types, unwrapped_coords, target_type):
    """Computes Mean Squared Dispacement curve (unit: ns, A^2)"""
    target_idx = np.nonzero(raw_types == target_type)[0]
    target_coords = unwrapped_coords[:, target_idx]

    ts = np.linspace(1, target_coords.shape[0] - 1, 100, dtype=int)
    msds = np.array([
        np.mean(np.sum((target_coords[t:] - target_coords[:-t])**2, axis=-1))
        for t in ts])
    # Convert to ns
    ts = ts * 2e-3
    return ts, msds


def get_conductivity(lattices, raw_types, unwrapped_coords, pop_mat):
    """Total conducitivty of the system (unit: S/cm)."""
    max_cluster = 10
    e_const = 1.6021766209e-19
    kb_const = 1.38064852e-23
    T = 353.0

    li_diff = get_diffusivity(
        raw_types, unwrapped_coords, target_type=90)  # cm^2/s
    tfsi_diff = get_diffusivity(
        raw_types, unwrapped_coords, target_type=93)  # cm^2/s

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


def get_cif_at_frame(wrapped_coords, lattices, atom_types, k):
    """Computes the cif text representation of the structure at frame k."""
    structure = Structure(
        lattice=np.diag(lattices[k]),
        species=atom_types,
        coords=wrapped_coords[k],
        coords_are_cartesian=True,
        to_unit_cell=True
    )
    return structure.to('cif')


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
