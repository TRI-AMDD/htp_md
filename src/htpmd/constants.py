"""
Module to store everything that has been constant for the past 13.8 billion years.
"""


ATOM_MASSES = [
    0.0, 1.008, 4.002602, 6.94, 9.0121831, 10.81, 12.011, 14.007, 15.999,
    18.998403163, 20.1797, 22.98976928, 24.305, 26.9815385, 28.085,
    30.973761998, 32.06, 35.45, 39.948, 39.0983, 40.078, 44.955908,
    47.867, 50.9415, 51.9961, 54.938044, 55.845, 58.933194, 58.6934, 63.546,
    65.38, 69.723, 72.63, 74.921595, 78.971, 79.904, 83.798, 85.4678, 87.62,
    88.90584, 91.224, 92.90637, 95.95, 97.90721, 101.07, 102.9055, 106.42,
    107.8682, 112.414, 114.818, 118.71, 121.76, 127.6, 126.90447, 131.293,
    132.90545196, 137.327, 138.90547, 140.116, 140.90766, 144.242, 144.91276,
    150.36, 151.964, 157.25, 158.92535, 162.5, 164.93033, 167.259, 168.93422,
    173.045, 174.9668, 178.49, 180.94788, 183.84, 186.207, 190.23, 192.217,
    195.084, 196.966569, 200.592, 204.38, 207.2, 208.9804, 209.0, 210.0,
    222.0, 223.0, 226.0, 227.0, 232.0377, 231.03588, 238.02891, 237.0, 244.0,
    243.0, 247.0, 247.0, 251.0, 252.0]


FARADAY_CONSTANT = 1.6021766209e-19
BOLTZMANN_CONSTANT = 1.38064852e-23
ATOMIC_MASS_IN_G = 1.66053906660e-24


ANGSTROM = 1e-10
CENTIMETER = 1e-2
NANOSECOND = 1e-9
PICOSECOND = 1e-12
KILOGRAM = 1e3


class AtomType:
    """Atomic numbers of different elements."""

    LITHIUM = 3
    NITROGEN = 7
    OXYGEN = 8
    SULFUR = 16


class RawType:
    """
    Raw atom type number used in LAMMPS.
    (Same element may have multiple raw types.)
    """

    LI = 90  # Li+
    TFSI = 93  # N atom in the TFSI-
