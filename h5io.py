import sys
import time

import numpy as np
import h5py
from htpmd.trajectory.load import LammpsTrajectoryLoader


def store_data(data_names, data, path):
    hf = h5py.File(path, 'w')
    for i in range(len(data_names)):
      hf.create_dataset(data_names[i], data=data[i])
    hf.close()


def load_data(data_names, path):
    hf = h5py.File(path, 'r')
    data = []
    for i in range(len(data_names)):
      d = np.array(hf.get(data_names[i]))
      data.append(d)
    hf.close()
    return data


def convert_to_h5(dir_name, out_file):
    t0 = time.time()
    traj = LammpsTrajectoryLoader().load(dir_name)
    t1 = time.time()
    print(f'Time for loading text files {t1 - t0}')
    # Just as an example, the list of data is incomplete here.
    store_data(
        ['lattices', 'unwrapped_coords', 'wrapped_coords'],
        [traj.lattices, traj.unwrapped_coords, traj.wrapped_coords],
        out_file,
    )


def load_h5(out_file):
    t0 = time.time()
    data = load_data(
        ['lattices', 'unwrapped_coords', 'wrapped_coords'],
        out_file
    )
    t1 = time.time()

    print(f'Time for loading h5 files {t1 - t0}')


def random_access(out_file):
    data = h5py.File(out_file)

    t0 = time.time()
    coord = data['unwrapped_coords'][2]
    t1 = time.time()

    print(f'Time for random access {t1 - t0}')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('[Usage] dir_name out_file')
        sys.exit(1)
    dir_name, out_file = sys.argv[1:]
    convert_to_h5(dir_name, out_file)
    load_h5(out_file)
    random_access(out_file)