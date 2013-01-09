import numpy as np

atomic_number = {
'H' : 1,                                                                                                                                                 'He': 2,
'Li': 3, 'Be': 4,                                                                                           'B' : 5, 'C' : 6, 'N' : 7, 'O' : 8, 'F' : 9, 'Ne':10,
'Na':11, 'Mg':12,                                                                                           'Al':13, 'Si':14, 'P' :15, 'S' :16, 'Cl':17, 'Ar':18,
'K' :19, 'Ca':20, 'Sc':21, 'Ti':22, 'V' :23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36,
'Rb':37, 'Sr':38, 'Y' :39, 'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49, 'Sn':50, 'Sb':51, 'Te':52, 'I' :53, 'Xe':54,
'Cs':55, 'Ba':56,          'Hf':72, 'Ta':73, 'W' :74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80, 'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86,
                }
class Geometry:
    def __init__(self, posfile=None, status='r'):
        self.posfile = posfile

        self.natoms = 0
        self.atom_types = list()  # Holds only the unique atoms
        self.atom_list = list()  # Holds all the atoms
        self.lattice_constant = 0.0
        self.unit_cell = np.zeros((3, 3))
        self.positions = None
        self.coordinate_type = None
        self.update_type = None
        self.number_atoms = 0
        self.element_numbers = None

        if (status == 'r'):
            self.read_posfile(posfile)

    def read_posfile(self, filename):
        """Reads in the geometry from a file"""
        raise NotImplementedError

    def write_posfile(self, filename):
        """Writes the geometry to a file"""
        raise NotImplementedError

#####################
# Utility functions #
#####################

    def direct2cart(self, positions):
        """Converts between direct coordinates and cartesian coordinates"""
        cell = self.unit_cell.copy() * self.lattice_constant
        for i in range(self.number_atoms):
            xnew = np.dot(cell[:, 0], positions[i])
            ynew = np.dot(cell[:, 1], positions[i])
            znew = np.dot(cell[:, 2], positions[i])

            positions[i] = np.array([xnew, ynew, znew])
        return positions

    def propogate_atoms(self, alist, nlist):
        """Given a list of unique atoms and a list of how many of each atom type
        is present, creates a list containing all the atoms"""
        typelist = list()

        for i in range(len(alist)):
            for j in range(nlist[i]):
                typelist.append(alist[i])

        return typelist

    def sort_atoms(self):
        """Sorts the atom list by atom type, and rearranges the positions
        match the new ordering"""
        sort_list = list()
        for i in range(len(self.atom_list)):
            atom = self.atom_list[i]
            anumber = atomic_number[atom]
            pos = self.positions[i]
            sort_list.append([atom, anumber, pos])

        new_positions = list()
        new_atom_list = list()
        sort_list = sorted(sort_list, key=lambda k: k[1])
        for sort in sort_list:
            new_atom_list.append(sort[0])
            new_positions.append(sort[2])

        self.set_positions(np.array(new_positions))
        self.set_atom_list(new_atom_list)


#################
# set functions #
#################
    def set_atom_types(self, atom_types):
        """Sets atom_types list"""
        self.atom_types = atom_types

    def set_unit_cell(self, cell):
        """Sets unit_cell np.array"""
        self.unit_cell = cell

    def set_positions(self, positions):
        """Sets positions np.array"""
        self.positions = positions

    def set_element_numbers(self, element_numbers):
        """Sets element_numbers list"""
        self.element_numbers = element_numbers

    def set_atom_list(self, atom_list):
        """Sets atom_list list"""
        self.atom_list = atom_list

#################
# get functions #
#################
    def get_atom_types(self):
        """Returns list of unique atom types"""
        return self.atom_types

    def get_unit_cell(self):
        """Returns the unit cell as a 3x3 np.array"""
        return self.unit_cell.copy()

    def get_positions(self):
        """Returns the atomic positions in cartesians as a Nx3 np.array"""
        return self.positions.copy()

    def get_element_numbers(self):
        """Returns a list of numbers for each element type """
        return self.element_numbers

    def get_atom_list(self):
        """Returns a list of each individual atom"""
        return self.atom_list

