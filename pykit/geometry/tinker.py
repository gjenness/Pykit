import os
import numpy as np

from pykit.geometry.geometry import Geometry

atomic_number = {
'H' : 1,                                                                                                                                                 'He': 2,
'Li': 3, 'Be': 4,                                                                                           'B' : 5, 'C' : 6, 'N' : 7, 'O' : 8, 'F' : 9, 'Ne':10,
'Na':11, 'Mg':12,                                                                                           'Al':13, 'Si':14, 'P' :15, 'S' :16, 'Cl':17, 'Ar':18,
'K' :19, 'Ca':20, 'Sc':21, 'Ti':22, 'V' :23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36,
'Rb':37, 'Sr':38, 'Y' :39, 'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49, 'Sn':50, 'Sb':51, 'Te':52, 'I' :53, 'Xe':54,
'Cs':55, 'Ba':56,          'Hf':72, 'Ta':73, 'W' :74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80, 'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86,
                }


class Tinker(Geometry):
    def __init__(self, posfile=None, status='r'):
        Geometry.__init__(self, posfile=posfile, status=status)

    def read_posfile(self, filename):
        try:
            posf = open(filename, 'r')
            contents = posf.readlines()
            posf.close()
        except:
            raise RuntimeError('File %s does not exist' % filename)

        self.natoms = int(contents[0].split()[0])

        alist = list()
        pos = list()

        for line in contents[1::]:
            data = line.split()
            alist.append(data[1])
            pos.append(data[2:5])

        self.atom_list = alist
        self.positions = np.array(pos, dtype='float')

    def write_posfile(self, filename):
        sym = self.get_atom_list()
        pos = self.get_positions().copy()

        posf = open(filename, 'w')
        posf.write('%i - Tinker geometry file prepared by PYKIT\n' % len(sym))
        for i in range(len(sym)):
            posf.write('%10i%10s' % (i+1, sym[i]))
            for j in range(3):
                posf.write('%16.8f' % pos[i, j])
            posf.write('%10i\n' % atomic_number[sym[i]])




