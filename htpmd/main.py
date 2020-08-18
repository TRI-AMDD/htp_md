"""

Usage:
  main.py <action> [-d <dir_path>]

Options:
  -h --help     Show this screen.
  -d --dir      Select directory of trajectory data

"""

from docopt import docopt
from htpmd import analysis
import pkg_resources


def analyze(dir_path):
    """
    Applies array of analyses on the trajectory data in the provided directory path and returns the results.

    Args:
        dir_path (str):                     path to trajectory data dictionary

    Returns:
        dict:                               analysis results

    """
    results = dict()
    # Load trajectories and pre-computed properties
    wrapped_coords, lattices, raw_types, atom_types, unwrapped_coords = (
        analysis.get_raw_traj(dir_path))
    pop_mat = analysis.get_population_matrix(dir_path)
    metadata = analysis.get_metadata(dir_path)

    # Compute properties
    results.update(metadata)
    results['li_diffusivity'] = analysis.get_diffusivity(
        raw_types, unwrapped_coords, target_type=90)
    results['tfsi_diffusivity'] = analysis.get_diffusivity(
        raw_types, unwrapped_coords, target_type=93))
    results['conductivity'] = analysis.get_conductivity(
        lattices, raw_types, unwrapped_coords, pop_mat)
    results['molarity'] =  analysis.get_molarity(raw_types, atom_types)
    return results


def get_version():
    """
    Returns the current htpmd module version.

    Returns:
        str:                htpmd version
    """

    version = pkg_resources.require("htpmd")[0].version
    return version


def main():

    arguments = docopt(__doc__)
    action = arguments['<action>']

    if action == 'analyze':

        if not arguments['--dir']:
            results = dict()
        else:
            dir_path = arguments['<dir_path>']
            results = analyze(dir_path)

    else:
        results = dict()

    results['htp_md_version'] = get_version()
    print(results)


if __name__ == '__main__':
    main()
