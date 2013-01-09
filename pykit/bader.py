import re
import os
import numpy as np

from pykit.geometry import geomread

nuclear_charges = {'H': 1, 'C': 4, 'O': 6, 'Rh': 9, 'Si': 4, 'Ti': 4, 'Pt': 10}


class bader:

# Use the initializer to set any files from a previous run
# run will run the Bader analysis
# get_summary will print a summary of results
# get_atom_charges will return an array of Bader charges

    def __init__(self, restart=False, chgcar='CHGCAR', chgcarsum='CHGCARsum',
                 datafile='ACF.dat'):
        self.restart = restart

        self.chgcar = chgcar
        self.chgcarsum = chgcarsum
        self.datafile = datafile

        self.atom_list = list()
        self.bader_charges = dict()

    def run(self):
#        chgsum = '/home/gjenness/programs/vtstscripts/chgsum.pl'
        baderprog = '/home/gjenness/programs/bader/bader'
        exitcode = os.system('%s %s -ref %s >& bader.out' % (baderprog, self.chgcar, self.chgcarsum))
        if (exitcode > 0):
            raise RuntimeError('Bader exited with error code %i' % exitcode)

    def read_bader_charges(self):
        charges = list()

        baderout = open(self.datafile, 'r')
        lines = baderout.readlines()

        for i in range(2, len(lines)-4):
            charges.append(float(lines[i].split()[4]))

        return np.array(charges)

    def calculate_atom_charges(self):
        znuc = list()
        baderchg = list()

        for i in range(len(self.atom_list)):
            znuc.append(nuclear_charges[self.atom_list[i]])

        baderchg = self.read_bader_charges()
        atom_charges = znuc - baderchg

        for i in range(len(self.atom_list)):
            label = self.atom_list[i] + str(i)
            self.bader_charges[label] = atom_charges[i]

    def position_sort(self, positions, coordinate='z'):
        col = {'x': 0, 'y': 1, 'z': 2}
        pos = positions[:, col[coordinate]]
        data = list()

        for i in range(len(self.atom_list)):
            atom = self.atom_list[i] + str(i)
            charge = self.bader_charges[atom]
            data.append([atom, pos[i], charge])

        data = sorted(data, key=lambda t: t[2])
        print data

    def set_atom_list(self, alist):
        self.atom_list = alist

    def get_bader_charges(self):
        self.read_bader_charges()
        return self.bader_charges

    def get_total_atom_charge(self, atomtype, exclude=None, include=None):
        chargesum = 0.0
        labellist = list()
        natoms = len(self.atom_list)

        if isinstance(atomtype, list):
            atomlist = atomtype
        else:
            atomlist = [atomtype]

        if isinstance(exclude, list):
            excludelist = exclude
        else:
            excludelist = [exclude]

        if isinstance(include, list):
            includelist = include
        else:
            includelist = [include]

        self.calculate_atom_charges()

        for atom in atomlist:
            for i in range(natoms):
                if (self.atom_list[i] == atom):
                    labellist.append(self.atom_list[i] + str(i))

        for key, val in self.bader_charges.items():
            if (key in labellist) and (key not in excludelist) or (key in includelist):
                chargesum += val

        return chargesum

    def get_summary(self):
        self.calculate_atom_charges()

        print 'Bader charge analysis results:\n'
        print '%-10s%10s%10s%20s' % ('Atom', 'ZNuc', 'Bader', 'Electron change?')
        print '--------------------------------------------------'
        for i in range(len(self.bader_charges)):
            atom = self.atom_list[i]
            znuc = nuclear_charges[atom]
            charge = self.bader_charges[atom + str(i)].copy()

            if (charge < 0.0):
                de = 'Gained electrons'
            elif (charge > 0.0):
                de = 'Lost electrons'
            elif (charge == 0.0):
                de = 'Remained the same'

            print '%-10s%10.4f%10.4f%20s' % (atom, znuc, charge, de)
        print '--------------------------------------------------'






