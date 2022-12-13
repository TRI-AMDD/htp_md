# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under
# the GNU General Public License.


oneline = "Read, write, manipulate dump files and particle attributes"

docstr = """
d = dump("dump.one")              read in one or more dump files
d = dump("dump.1 dump.2.gz")	  can be gzipped
d = dump("dump.*")		  wildcard expands to multiple files
d = dump("dump.*",0)		  two args = store filenames, but don't read

  incomplete and duplicate snapshots are deleted
  atoms will be unscaled if stored in files as scaled
  self-describing column names assigned 

time = d.next()             	  read next snapshot from dump files

  used with 2-argument constructor to allow reading snapshots one-at-a-time
  snapshot will be skipped only if another snapshot has same time stamp
  return time stamp of snapshot read
  return -1 if no snapshots left or last snapshot is incomplete
  no column name assignment or unscaling is performed

d.map(1,"id",3,"x")               assign names to columns (1-N)

  not needed if dump file is self-describing

d.tselect.all()			  select all timesteps
d.tselect.one(N)		  select only timestep N
d.tselect.none()		  deselect all timesteps
d.tselect.skip(M)		  select every Mth step
d.tselect.test("$t >= 100 and $t < 10000")      select matching timesteps
d.delete()	      	      	  delete non-selected timesteps

  selecting a timestep also selects all atoms in the timestep
  skip() and test() only select from currently selected timesteps
  test() uses a Python Boolean expression with $t for timestep value
    Python comparison syntax: == != < > <= >= and or

d.aselect.all()	      	                      select all atoms in all steps
d.aselect.all(N)      	                      select all atoms in one step
d.aselect.test("$id > 100 and $type == 2")    select match atoms in all steps
d.aselect.test("$id > 100 and $type == 2",N)  select matching atoms in one step

  all() with no args selects atoms from currently selected timesteps
  test() with one arg selects atoms from currently selected timesteps
  test() sub-selects from currently selected atoms
  test() uses a Python Boolean expression with $ for atom attributes
    Python comparison syntax: == != < > <= >= and or
    $name must end with a space

d.write("file")	   	           write selected steps/atoms to dump file
d.write("file",head,app)	   write selected steps/atoms to dump file
d.scatter("tmp")		   write selected steps/atoms to multiple files

  write() can be specified with 2 additional flags
    head = 0/1 for no/yes snapshot header, app = 0/1 for write vs append
  scatter() files are given timestep suffix: e.g. tmp.0, tmp.100, etc

d.scale() 	    	  	   scale x,y,z to 0-1 for all timesteps
d.scale(100)			   scale atom coords for timestep N
d.unscale()			   unscale x,y,z to box size to all timesteps
d.unscale(1000)			   unscale atom coords for timestep N
d.wrap()			   wrap x,y,z into periodic box via ix,iy,iz
d.unwrap()			   unwrap x,y,z out of box via ix,iy,iz
d.owrap("other")		   wrap x,y,z to same image as another atom
d.sort()              	  	   sort atoms by atom ID in all selected steps
d.sort("x")            	  	   sort atoms by column value in all steps
d.sort(1000)			   sort atoms in timestep N
t = d.time()  	     	       	   return vector of selected timestep values
"""

import sys, subprocess, re, glob, types
from os import popen
from math import *  # any function could be used by set()

try:
    import numpy as np

    oldnumeric = False
except:
    import Numeric as np

    oldnumeric = True


# try: from DEFAULTS import PIZZA_GUNZIP
# except: PIZZA_GUNZIP = "gunzip"

# Class definition

class dump:

    # --------------------------------------------------------------------

    def __init__(self, *list):
        self.snaps = []
        self.nsnaps = self.nselect = 0
        self.names = {}
        self.tselect = tselect(self)
        self.aselect = aselect(self)
        self.atype = "type"
        self.bondflag = 0
        self.bondlist = []
        self.triflag = 0
        self.trilist = []
        self.lineflag = 0
        self.linelist = []
        self.objextra = None

        # flist = list of all dump file names

        words = list[0].split()
        self.flist = []
        for word in words: self.flist += glob.glob(word)
        if len(self.flist) == 0 and len(list) == 1:
            raise Exception("no dump file specified")

        if len(list) == 1:
            self.increment = 0
            self.read_all()
        else:
            self.increment = 1
            self.nextfile = 0
            self.eof = 0

    # --------------------------------------------------------------------

    def read_all(self):

        # read all snapshots from each file
        # test for gzipped files

        for file in self.flist:
            if file[-3:] == ".gz":
                f = popen("%s -c %s" % (PIZZA_GUNZIP, file), 'r')
            else:
                f = open(file)

            snap = self.read_snapshot(f)
            while snap:
                self.snaps.append(snap)
                print(snap.time, end=' ')
                sys.stdout.flush()
                snap = self.read_snapshot(f)

            f.close()
        print()

        # sort entries by timestep, cull duplicates

        # self.snaps.sort(self.compare_time)
        self.cull()
        self.nsnaps = len(self.snaps)
        print("read %d snapshots" % self.nsnaps)

        # select all timesteps and atoms

        self.tselect.all()

        # print column assignments

        if len(self.names):
            print("assigned columns:", self.names2str())
        else:
            print("no column assignments made")

        # if snapshots are scaled, unscale them

        if ("x" not in self.names) or \
                ("y" not in self.names) or \
                ("z" not in self.names):
            print("dump scaling status is unknown")
        elif self.nsnaps > 0:
            if self.scale_original == 1:
                self.unscale()
            elif self.scale_original == 0:
                print("dump is already unscaled")
            else:
                print("dump scaling status is unknown")

    def read_snapshot(self, f):
        try:
            snap = Snap()
            item = f.readline()
            snap.time = int(f.readline().split()[0])  # just grab 1st field
            item = f.readline()
            snap.natoms = int(f.readline())

            snap.aselect = np.zeros(snap.natoms)

            item = f.readline()
            words = item.split("BOUNDS ")
            if len(words) == 1:
                snap.boxstr = ""
            else:
                snap.boxstr = words[1].strip()
            if "xy" in snap.boxstr:
                snap.triclinic = 1
            else:
                snap.triclinic = 0

            words = f.readline().split()
            if len(words) == 2:
                snap.xlo, snap.xhi, snap.xy = float(words[0]), float(words[1]), 0.0
            else:
                snap.xlo, snap.xhi, snap.xy = \
                    float(words[0]), float(words[1]), float(words[2])

            words = f.readline().split()
            if len(words) == 2:
                snap.ylo, snap.yhi, snap.xz = float(words[0]), float(words[1]), 0.0
            else:
                snap.ylo, snap.yhi, snap.xz = \
                    float(words[0]), float(words[1]), float(words[2])

            words = f.readline().split()
            if len(words) == 2:
                snap.zlo, snap.zhi, snap.yz = float(words[0]), float(words[1]), 0.0
            else:
                snap.zlo, snap.zhi, snap.yz = \
                    float(words[0]), float(words[1]), float(words[2])

            item = f.readline()
            if len(self.names) == 0:
                self.scale_original = -1
                xflag = yflag = zflag = -1
                words = item.split()[2:]
                if len(words):
                    for i in range(len(words)):
                        if words[i] == "x" or words[i] == "xu":
                            xflag = 0
                            self.names["x"] = i
                        elif words[i] == "xs" or words[i] == "xsu":
                            xflag = 1
                            self.names["x"] = i
                        elif words[i] == "y" or words[i] == "yu":
                            yflag = 0
                            self.names["y"] = i
                        elif words[i] == "ys" or words[i] == "ysu":
                            yflag = 1
                            self.names["y"] = i
                        elif words[i] == "z" or words[i] == "zu":
                            zflag = 0
                            self.names["z"] = i
                        elif words[i] == "zs" or words[i] == "zsu":
                            zflag = 1
                            self.names["z"] = i
                        else:
                            self.names[words[i]] = i
                    if xflag == 0 and yflag == 0 and zflag == 0:
                        self.scale_original = 0
                    if xflag == 1 and yflag == 1 and zflag == 1:
                        self.scale_original = 1

            if snap.natoms:
                words = f.readline().split()
                ncol = len(words)
                for i in range(1, snap.natoms):
                    words += f.readline().split()
                floats = list(map(float, words))
                if oldnumeric:
                    atoms = np.zeros((snap.natoms, ncol), np.Float)
                else:
                    atoms = np.zeros((snap.natoms, ncol), np.float)
                start = 0
                stop = ncol
                for i in range(snap.natoms):
                    atoms[i] = floats[start:stop]
                    start = stop
                    stop += ncol
            else:
                atoms = None
            snap.atoms = atoms
            return snap
        except:
            return 0

    # --------------------------------------------------------------------
    # unscale coords from 0-1 to box size for all snapshots or just one
    # use 6 params as h-matrix to treat orthogonal or triclinic boxes

    def unscale(self, *list):
        if len(list) == 0:
            print("Unscaling dump ...")
            x = self.names["x"]
            y = self.names["y"]
            z = self.names["z"]
            for snap in self.snaps:
                self.unscale_one(snap, x, y, z)
        else:
            i = self.findtime(list[0])
            x = self.names["x"]
            y = self.names["y"]
            z = self.names["z"]
            self.unscale_one(self.snaps[i], x, y, z)

    # --------------------------------------------------------------------

    def unscale_one(self, snap, x, y, z):
        if snap.xy == 0.0 and snap.xz == 0.0 and snap.yz == 0.0:
            xprd = snap.xhi - snap.xlo
            yprd = snap.yhi - snap.ylo
            zprd = snap.zhi - snap.zlo
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = snap.xlo + atoms[:, x] * xprd
                atoms[:, y] = snap.ylo + atoms[:, y] * yprd
                atoms[:, z] = snap.zlo + atoms[:, z] * zprd
        else:
            xlo_bound = snap.xlo
            xhi_bound = snap.xhi
            ylo_bound = snap.ylo
            yhi_bound = snap.yhi
            zlo_bound = snap.zlo
            zhi_bound = snap.zhi
            xy = snap.xy
            xz = snap.xz
            yz = snap.yz
            xlo = xlo_bound - min((0.0, xy, xz, xy + xz))
            xhi = xhi_bound - max((0.0, xy, xz, xy + xz))
            ylo = ylo_bound - min((0.0, yz))
            yhi = yhi_bound - max((0.0, yz))
            zlo = zlo_bound
            zhi = zhi_bound
            h0 = xhi - xlo
            h1 = yhi - ylo
            h2 = zhi - zlo
            h3 = yz
            h4 = xz
            h5 = xy
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = snap.xlo + atoms[:, x] * h0 + atoms[:, y] * h5 + atoms[:, z] * h4
                atoms[:, y] = snap.ylo + atoms[:, y] * h1 + atoms[:, z] * h3
                atoms[:, z] = snap.zlo + atoms[:, z] * h2

    # --------------------------------------------------------------------
    # convert column names assignment to a string, in column order

    def names2str(self):
        pairs = list(self.names.items())
        values = list(self.names.values())
        ncol = len(pairs)
        str = ""
        for i in range(ncol):
            if i in values:
                str += pairs[values.index(i)][0] + ' '
        return str

    # --------------------------------------------------------------------
    # sort atoms by atom ID in all selected timesteps by default
    # if arg = string, sort all steps by that column
    # if arg = numeric, sort atoms in single step

    def sort(self, *list):
        if len(list) == 0:
            print("Sorting selected snapshots ...")
            id = self.names["id"]
            for snap in self.snaps:
                if snap.tselect:
                    self.sort_one(snap, id)
        elif type(list[0]) is bytes:
            print("Sorting selected snapshots by %s ..." % list[0])
            id = self.names[list[0]]
            for snap in self.snaps:
                if snap.tselect:
                    self.sort_one(snap, id)
        else:
            i = self.findtime(list[0])
            id = self.names["id"]
            self.sort_one(self.snaps[i], id)

    # --------------------------------------------------------------------
    # sort a single snapshot by ID column

    def sort_one(self, snap, id):
        atoms = snap.atoms
        ids = atoms[:, id]
        ordering = np.argsort(ids)
        for i in range(len(atoms[0])):
            atoms[:, i] = np.take(atoms[:, i], ordering)

    # --------------------------------------------------------------------
    # return vector of selected snapshot time stamps

    def time(self):
        vec = self.nselect * [0]
        i = 0
        for snap in self.snaps:
            if not snap.tselect:
                continue
            vec[i] = snap.time
            i += 1
        return vec

    # --------------------------------------------------------------------
    # extract vector(s) of values for selected atoms at chosen timestep

    def vecs(self, n, *list):
        snap = self.snaps[self.findtime(n)]

        if len(list) == 0:
            raise Exception("no columns specified")
        columns = []
        values = []
        for name in list:
            columns.append(self.names[name])
            values.append(snap.nselect * [0])
        ncol = len(columns)

        m = 0
        for i in range(snap.natoms):
            if not snap.aselect[i]:
                continue
            for j in range(ncol):
                values[j][m] = snap.atoms[i][columns[j]]
            m += 1

        if len(list) == 1:
            return values[0]
        else:
            return values

    # --------------------------------------------------------------------
    # sort snapshots on time stamp

    def compare_time(self, a, b):
        if a.time < b.time:
            return -1
        elif a.time > b.time:
            return 1
        else:
            return 0

    # --------------------------------------------------------------------
    # delete successive snapshots with duplicate time stamp

    def cull(self):
        i = 1
        while i < len(self.snaps):
            if self.snaps[i].time == self.snaps[i - 1].time:
                del self.snaps[i]
            else:
                i += 1

    #     return time,box,atoms,bonds,tris,lines

    # --------------------------------------------------------------------

    def findtime(self, n):
        for i in range(self.nsnaps):
            if self.snaps[i].time == n:
                return i
        raise Exception("no step %d exists" % n)


# --------------------------------------------------------------------
# one snapshot

class Snap:
    pass


# --------------------------------------------------------------------
# time selection class

class tselect:

    def __init__(self, data):
        self.data = data

    # --------------------------------------------------------------------

    def all(self):
        data = self.data
        for snap in data.snaps:
            snap.tselect = 1
        data.nselect = len(data.snaps)
        data.aselect.all()
        print("%d snapshots selected out of %d" % (data.nselect, data.nsnaps))

    # --------------------------------------------------------------------

    def one(self, n):
        data = self.data
        for snap in data.snaps:
            snap.tselect = 0
        i = data.findtime(n)
        data.snaps[i].tselect = 1
        data.nselect = 1
        data.aselect.all()
        print("%d snapshots selected out of %d" % (data.nselect, data.nsnaps))

    # --------------------------------------------------------------------

    def none(self):
        data = self.data
        for snap in data.snaps:
            snap.tselect = 0
        data.nselect = 0
        print("%d snapshots selected out of %d" % (data.nselect, data.nsnaps))


# # --------------------------------------------------------------------
# # atom selection class

class aselect:

    def __init__(self, data):
        self.data = data

    # --------------------------------------------------------------------

    def all(self, *args):
        data = self.data
        if len(args) == 0:  # all selected timesteps
            for snap in data.snaps:
                if not snap.tselect:
                    continue
                for i in range(snap.natoms):
                    snap.aselect[i] = 1
                snap.nselect = snap.natoms
        else:  # one timestep
            n = data.findtime(args[0])
            snap = data.snaps[n]
            for i in range(snap.natoms):
                snap.aselect[i] = 1
            snap.nselect = snap.natoms
