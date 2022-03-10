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
    return analysis.get_all_properties(dir_path)


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
            raise SyntaxError("Directory Not Provided")
            results = dict()
        else:
            dir_path = arguments['<dir_path>']
            results = analyze(dir_path)

    else:
        raise SyntaxError("Invalid Action")
        results = dict()

    results['htp_md_version'] = get_version()
    print(results)


if __name__ == '__main__':
    main()


