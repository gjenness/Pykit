import os, math
import numpy as np

from pykit.geometry import geomread
from pykit.utility import convert2list, integrate


class DensityOfStates:
    """
    dosfile is the file containing the DOS data
    posfile is the file containing the geometry

    energies and DOS are stored in a dictionary called self.dos_data
    with the energies and DOS being accessed via one of two keys:
        'all':  returns the DOS of the whole system
        atom label + atom index:  returns the DOS for that particular atom
    when dos_data is accessed, it'll return a list with two entries.  The first entry
    will contain the energies, and the 2nd entry will contain the DOS
    """

    def __init__(self, dosfile=None, posfile=None, status=None):
#        print 'DensityOfStates.__init__() called'
        self.natoms = 0    # Number of atoms
        self.efermi = 0.0  # Fermi energy
        self.npoints = 0   # Number of points print for the DOS
        self.spin = False  # Is spin on?
        self.pdos = False  # PDOS?
        self.opdos = False # Orbital PDOS?
        self.energies = list()  # Holds the energy grid
        self.dos_data = dict()  # Holds all the information regarding the DOS

        self.atom_types = list()  # Stores the atom types in the dosfile, taken from a geometry file
        self.atom_labels = list() # Stores the atom types with an integer label

#######################################################
# Next two functions are defined by the derived class #
#######################################################
    def read_dosfile(self, dosfile):
        raise NotImplementedError

    def process_dos_data(self, lines):
        raise NotImplementedError

####################################################################################
# Remaining class functions should not be changed when creating the derived class! #
####################################################################################

# calculate_band_energy returns a weighted average of the band energy for a given DOS
# by default, analyzes the d-band for the whole system
    def calculate_band_energy(self, band='d', atomtype=None, atomindex=None,
                              exclude=None, include=None):
        if (atomtype is not None) and (atomindex is not None):
            raise RuntimeError('atomtype and atomindex are mutually exclusive!')

        if (atomtype is None) and (atomindex is None):
            print 'Analyzing the %s-band for the whole system' % band.lower()
            atomtype = 'all'

        if (atomtype is not 'all') and (not self.pdos):
            raise RuntimeError('Cannot analyze %s atom(s), PDOS data not found.' % atomtype)

        atomlist = convert2list(atomtype)
        excludelist = convert2list(exclude)
        includelist = convert2list(include)

        allowed_bands = ['s', 'p', 'd']  # List of allowed bands
        nbands = len(band)  # nbands is the number of bands to be analyzed
        bandlist = list()  # List of bands to be analyzed

        for ban in band:
            if (ban.lower() not in allowed_bands):
                raise RuntimeError('Unknown band type: %s' % ban)
            else:
                bandlist.append(ban.lower())

        statesum = 0.0  # stores the sum of the states
        esum = 0.0  # stores the sum of the energy weighted by the number of states

        if (atomtype == 'all'):
            energy, dos = self.get_total_dos()
            esum, statesum = self.calculate_weighted_average(energy, dos, band)
        else:
            if (not atomlist) and (not includelist):
                raise RuntimeError('No atoms specified!')
            elif (not atomlist) and (includelist):
                dos = self.get_total_atom_dos(atomtype=None, exclude=excludelist,
                                                      include=includelist)
                es, ss = self.calculate_weighted_average(dos, band)
                esum += es
                statesum += ss
            else:
                for atom in atomlist:
                    dos = self.get_total_atom_dos(atomtype=atom, exclude=excludelist,
                                                  include=includelist)
                    es, ss = self.calculate_weighted_average(dos, band)
                    esum += es
                    statesum += ss

#        print 'Total:', esum, statesum
        return esum / statesum

    def calculate_weighted_average(self, dos=None, band=None):
        labels = {'s': 0, 'p': 1, 'd': 2}
        spinlabels = {'s': [0, 1], 'p': [2, 3], 'd': [4, 5]}
        energy = self.energies.copy()  # Use a copy of the energy grid just to be safe

        statesum = 0.0 # stores the sum of the states
        esum = 0.0  # stores the sum of the energy weighted by the number of states
        nbands = len(band)

        if (self.spin):
            for b in band:
                colup = spinlabels[b][0]
                coldown = spinlabels[b][1]
                for j in range(self.npoints):
                    statesum += (dos[j, colup] + dos[j, coldown])
                    esum += (self.energies[j] * (dos[j, colup] + dos[j, coldown]))
        elif (not self.spin):
            for b in band:
                col = labels[b]
                eweight = energy * dos[:, col]
                statesum = sum(dos[:, col])
                esum = sum(eweight)

        return (esum, statesum)

#####################
# Utility Functions #
#####################
    def help(self):
        print '\n\npykit.densityofstates help information:'
        print '\nset functions:'
        print '    set_atom_types(self, atomtypes): sets atom types'
        print '\nget functions:'
        print '    get_total_dos(self): returns energies and dos for whole system'
        print '    get_atom_dos(self, label): return energies and dos for a single atom'
        print '    get_total_atom_dos(self, atomtype): return energies and dos for a single atom type'
        print '    get_band_energy(self, band, atomtype, atomindex): returns a weighted averge band energy'
        print '\n\n'

    def reset_fermi_level(self, fermi):
        print 'Re-setting energy grid to new Fermi level.'
        self.energies += (self.efermi - fermi)
        print 'Fermi energy is now', fermi, 'eV.'
        self.__set_fermi_level(fermi)

#################
# Set Functions #
#################
    def set_atom_types(self, atomtypes):
        self.atom_types = atomtypes
        for i in range(len(self.atom_types)):
            self.atom_labels.append(self.atom_types[i] + str(i))

    def __set_fermi_level(self, fermi):
        """
        Only use if you know what you're doing!
        """
        self.efermi = fermi

#################
# Get Functions #
#################
    def get_energy_grid(self):
        """Returns the energy grid the DOS was calculated on"""
        return self.energies.copy()

    def get_total_dos(self):
        energies = self.dos_data['total'][0].copy()
        dos = self.dos_data['total'][1].copy()

        return (energies, dos)

    def get_atom_dos(self, label):
        """Given an atom label, returns the DOS for 
        that particular atom"""
        return self.dos_data[label][1].copy()

    def get_total_atom_dos(self, atomtype=None, exclude=list(), include=list()):
        """Given a particular atom type, returns the energies and total DOS 
        for that atom"""
        if (self.spin):
            totaldos = np.zeros((self.npoints, 6)) 
        elif (not self.spin):
            totaldos = np.zeros((self.npoints, 3)) 

        if (atomtype == None) or (not atomtype):
            atomtype = ' ' #  Turn atomtype into a blank string

        for label in self.atom_labels:
            if ((atomtype in label) or (label in include)) and (label not in exclude):
                totaldos += self.get_atom_dos(label).copy()

        return totaldos

    def get_band_energy(self, band='d', atomtype=None, atomindex=None,
                        exclude=None, include=None):
        return self.calculate_band_energy(band, atomtype=atomtype, atomindex=atomindex,
                                          exclude=exclude, include=include)

# Need to clean up the next to class functions and remove the duplicate code
# But for right now, we'll leave it as this, 'cause I'm lazy.
    def get_band_filling(self, band='d', atomtype=None, exclude=None,
                         include=None, erange=None):
        labels = {'s': 0, 'p': 1, 'd': 2}
        spinlabels = {'s': [0, 1], 'p': [2, 3], 'd': [4, 5]}
        allowed_bands = ['s', 'p', 'd']  # List of allowed bands
        bandlist = list()  # List of bands to be analyzed

        if (self.spin):
            totaldos = np.zeros((self.npoints, 6))
        elif (not self.spin):
            totaldos = np.zeros((self.npoints, 3)) 

        dosband = np.zeros((self.npoints))

        atomlist = convert2list(atomtype)
        excludelist = convert2list(exclude)
        includelist = convert2list(include)

        for ban in band:
            if (ban.lower() not in allowed_bands):
                raise RuntimeError('Unknown band type: %s' % ban)
            else:
                bandlist.append(ban.lower())

# Here we get the totaldos for each atom in the system
        for label in self.atom_labels:
            for atom in atomlist:
                if ((atom in label) or (label in includelist)) and (label not in excludelist):
                    totaldos += self.get_total_atom_dos(atomtype=atom,
                                                        exclude=excludelist,
                                                        include=includelist)

        if (erange == None):
#            print 'Integrating over whole range up to the Fermi level.'
            emin = self.energies[0]
            emax = 0.0
        elif (erange[1] == 'fermi'):
            emin = erange[0]
            emax = 0.0
        else:
            emin = erange[0]
            emax = erange[1]

        energylist = list()
        indexlist = list()

        for n,e in enumerate(self.energies):
            if (e > emin) and (e < emax):
                energylist.append(e)
                indexlist.append(n)

        start = indexlist[0]
        finish = indexlist[-1]
        if (self.spin):
            for b in band:
                colup = spinlabels[b][0]
                coldown = spinlabels[b][1]
                dosband += (totaldos[:, colup] + totaldos[:, coldown]) 
        elif (not self.spin):
            for b in band:
                col = labels[b]
                dosband += totaldos[:, col]

        trimdos = dosband[start:finish+1]

        return integrate(grid=energylist, data=trimdos)
        
    def get_band_width(self, band='d', atomtype=None, exclude=None,
                         include=None, erange=None):
# Yes I'm lazy, sue me.
        labels = {'s': 0, 'p': 1, 'd': 2}
        spinlabels = {'s': [0, 1], 'p': [2, 3], 'd': [4, 5]}
        allowed_bands = ['s', 'p', 'd']  # List of allowed bands
        bandlist = list()  # List of bands to be analyzed

        if (self.spin):
            totaldos = np.zeros((self.npoints, 6))
        elif (not self.spin):
            totaldos = np.zeros((self.npoints, 3)) 

        dosband = np.zeros((self.npoints))

        atomlist = convert2list(atomtype)
        excludelist = convert2list(exclude)
        includelist = convert2list(include)

        for ban in band:
            if (ban.lower() not in allowed_bands):
                raise RuntimeError('Unknown band type: %s' % ban)
            else:
                bandlist.append(ban.lower())

# Here we get the totaldos for each atom in the system
        for label in self.atom_labels:
            for atom in atomlist:
                if ((atom in label) or (label in includelist)) and (label not in excludelist):
                    totaldos += self.get_total_atom_dos(atomtype=atom,
                                                        exclude=excludelist,
                                                        include=includelist)

        if (erange == None):
#            print 'Integrating over whole range up to the Fermi level.'
            emin = self.energies[0]
            emax = 0.0
        elif (erange[1] == 'fermi'):
            emin = erange[0]
            emax = 0.0
        else:
            emin = erange[0]
            emax = erange[1]

        energylist = list()
        indexlist = list()

        for n,e in enumerate(self.energies):
            if (e > emin) and (e < emax):
                energylist.append(e)
                indexlist.append(n)

        start = indexlist[0]
        finish = indexlist[-1]
        if (self.spin):
            for b in band:
                colup = spinlabels[b][0]
                coldown = spinlabels[b][1]
                dosband += (totaldos[:, colup] + totaldos[:, coldown]) 
        elif (not self.spin):
            for b in band:
                col = labels[b]
                dosband += totaldos[:, col]

        trimdos = dosband[start:finish+1] * energylist * energylist

        filling = self.get_band_filling(band=band, atomtype=atomtype,
                                      include=include, exclude=exclude, erange=erange)

        return math.sqrt(integrate(grid=energylist, data=trimdos) / filling)
