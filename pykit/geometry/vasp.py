import os, itertools
import numpy as np

from pykit.geometry.geometry import Geometry


class Poscar(Geometry):
    def __init__(self, posfile=None, status='r'):
        Geometry.__init__(self, posfile=posfile, status=status)

    def read_posfile(self, filename):
# Just assume VASP4 format, and the title has the atom types as an ASE list
        try:
            posf = open(filename, 'r')
            contents = posf.readlines()
            posf.close()
        except:
            raise RuntimeError('File %s does not exist' % filename)

        self.atom_types = contents[0].split()
        self.lattice_constant = float(contents[1])

        cell = list()
        for i in range(3):
            line = contents[2 + i].split()
            tmp = [float(line[0]), float(line[1]), float(line[2])]
            cell.append(tmp)
        self.unit_cell = np.array(cell)

        self.element_numbers = contents[5].split()

        line6 = contents[6].strip().lower()
        if ('direct' in line6) or ('cart' in line6):
            self.coordinate_type = line6
            geomstart = 7
        else:
            self.update_type = line6
            self.coordinate_type = contents[7].strip().lower()
            geomstart = 8

#        print self.coordinate_type
        for i in range(len(self.element_numbers)):
            self.element_numbers[i] = int(self.element_numbers[i])

        self.number_atoms = sum(self.element_numbers)
        self.atom_list = self.propogate_atoms(self.atom_types, self.element_numbers)

        self.positions = np.zeros((self.number_atoms, 3))
        for i in range(self.number_atoms):
            tmp = contents[geomstart + i].split()
            for j in range(3):
                self.positions[i, j] = float(tmp[j])

        if (self.coordinate_type == 'direct'):
            self.positions = self.direct2cart(self.positions.copy())

    def write_posfile(self, filename):
        self.sort_atoms()

        unique_atoms = list()
        for list_, count in itertools.groupby(self.atom_list):
            number = 0
            for c in count:
                number += 1
            unique_atoms.append([list_, number])

        fpos = open(filename, 'w')
        for atom in unique_atoms:
            fpos.write('%s ' % atom[0])

        if (self.lattice_constant == 0.0):
            fpos.write('\n1.00\n')
        else:
            fpos.write('\n%f\n' % self.lattice_constant)

        for vec in self.unit_cell:
            for i in range(3):
                fpos.write('%18.8f' % vec[i])
            fpos.write('\n')

        for atom in unique_atoms:
            fpos.write('%10i' % atom[1])

        fpos.write('\ncartesian\n')

        for pos in self.positions:
            for i in range(3):
                fpos.write('%18.8f' % pos[i])
            fpos.write('\n')
        fpos.write('\n')
