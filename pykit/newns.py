import numpy as np
import numpy.linalg as la
import math


allowed_keywords = ['adsorbate_energies',
                    'adsorbate_electrons',
                    'continuum_range',
                    'kstates',
                    'continuum_electrons',
                    'metal_adsorbate_coupling',
                    'U_repulsion',
                    'scf_cycles',
                    'scf_thresh',
                   ]

class Newns:
    def __init__(self, **kwargs):
        print 'Solver for the Newns-Anderson model of chemisorption.'

# dict containing all the information about the orbitals
        self.orbital_information = dict()

# dict containing all the keywords and their values.
        self.keywords = dict()
# Set the keywords based on what is given in **kwargs
        for keyword in allowed_keywords:
            self.keywords[keyword] = None
        self.__set_keywords(**kwargs)

#        print self.keywords
# If true, self.run() needs to be given before anything else can be done
        self.need_run = True

    def __set_keywords(self, **kwargs):
        """Sets the keywords from kwargs"""
        for keyword in kwargs:
            if self.keywords.has_key(keyword):
                self.keywords[keyword] = kwargs[keyword]
 
        print '\nThe following input options were given:'
        for key, val in self.keywords.items():
            if (val is not None):
                print '%4s%-25s:%15s' % (' ', str(key), str(val))

        print 'Setting appropiate defaults for non-specified keywords.'
        self.__set_defaults()

        print '\nUse Newns.run() to generate the required Hamiltonian and orbital energies.'

    def __set_defaults(self):
        """Sets some default values"""
        if (self.keywords['scf_cycles'] == None):
            self.keywords['scf_cycles'] = 30
        if (self.keywords['U_repulsion'] == None):
            self.keywords['U_repulsion'] = 0.0
        if (self.keywords['scf_thresh'] == None):
            self.keywords['scf_thresh'] = 1e-2

    def __assign_electrons(self, total_electrons):
        """Assigns electrons according to alpha and beta spin.  By convention,
           I assume that the odd electron belongs to alpha spin"""
        if ((total_electrons % 2.0) > 0.0):
            alpha_electrons = int(math.floor(total_electrons / 2))
            beta_electrons = alpha_electrons
            alpha_electrons += (total_electrons - (alpha_electrons + beta_electrons))
        elif ((total_electrons % 2.0) == 0.0):
            alpha_electrons = int(math.floor(total_electrons / 2))
            beta_electrons = alpha_electrons
        elif (total_electrons == 1):
            alpha_electrons = 1
            beta_electrons = 0

        return (alpha_electrons, beta_electrons)

    def __build_Hamiltonian(self, ads_energy, ads_ele, kstates, con_ele, vak):
        """Builds the Newns-Anderson Hamiltonian, neglecting the e-e repulsion
           term.  This term is added back in when we consider a specific spin
           Hamiltonian"""

        las = len(ads_energy)  # Number of adsorbate orbitals
        lss = len(kstates)  # Number of continuum states

        nelements = las + lss
        hamiltonian = np.zeros( (nelements, nelements) )

# For debugging purposes
#        print len(sstates)
#        print sstates

# Rock that diagonal!
# Some counters
        nas = 0  # Counts the number of adsorbate orbitals
        nss = 0  # Counts the number of continuum states

        for i in range(nelements):
            if (i < las):
                hamiltonian[i, i] = ads_energy[nas]
                nas += 1
            elif (i < (las + lss)):
                hamiltonian[i, i] = kstates[nss]
                nss += 1

        # Rock that off-diagonal!
        for i in range(las, nelements):
            for j in range(las):
                if (i != j):
                    hamiltonian[i, j] = vak

                hamiltonian[j, i] = hamiltonian[i, j]

# DEBUG:  writes the hamiltonian to a text file
#        txt = open('hamil.txt', 'w')
#        for i in range(nelements):
#            for j in range(nelements):
#                txt.write('%9.3f' % hamiltonian[i, j])
#            txt.write('\n')
#        txt.close()

        return hamiltonian

# Keep an eye on this function!
    def __update_hamiltonian(self, Uold, Unew, hamiltonian):
        """Updates the Hamiltonian"""
        nads = len(self.keywords['adsorbate_energies'])
        diff = Unew - Uold
        for i in range(nads):
            hamiltonian[i, i] += diff

        return hamiltonian

    def __update_repulsion(self, Urep, coeff, electrons):
        """Sets the adsorbate e-e repulsion"""
        occupancy = 0.0
#        print coeff
        for i in range(electrons):
#            print coeff[0, i]
            occupancy += coeff[0, i]**2

        return (Urep*occupancy, occupancy)

    def __sort_orbitals(self, orbitals, orbital_coeff):
        """Sorts the orbitals and coefficients"""
        index = orbitals.argsort()

        orbitals = orbitals[index]
        orbital_coeff = orbital_coeff[:, index]

#        print '\nAfter sort:'
#        print orbitals
#        print orbital_coeff

        return (orbitals, orbital_coeff)
    
    def __scf(self, alpha_hamiltonian, beta_hamiltonian):
        """Performs the SCF"""
        print '\nPerforming SCF iterations.'
        print 'Grabbing relevant orbital information.'

        alpha_electrons = int(self.orbital_information['alpha_electrons'])
        beta_electrons = int(self.orbital_information['beta_electrons'])
        Uold_alpha = float(self.keywords['U_repulsion'])
        Uold_beta = float(self.keywords['U_repulsion'])

        ncycles = int(self.keywords['scf_cycles'])
        thresh = float(self.keywords['scf_thresh'])
        print '%-20s:%10i' % ('Max SCF cycles', ncycles)
        print '%-20s:%10.1e\n' % ('SCF Threshold', thresh)

        for i in range(ncycles):
            print 'Iteration %-5i' % i
            alpha_orbitals, alpha_coeff = la.eig(alpha_hamiltonian)
            beta_orbitals, beta_coeff = la.eig(beta_hamiltonian)

#            print 'Sorting alpha:'
            alpha_orbitals, alpha_coeff = self.__sort_orbitals(alpha_orbitals, alpha_coeff)
#            print 'Sorting beta:'
            beta_orbitals, beta_coeff = self.__sort_orbitals(beta_orbitals, beta_coeff)

            Unew_alpha, occ_alpha = self.__update_repulsion(Uold_alpha, beta_coeff , beta_electrons)
            Unew_beta , occ_beta  = self.__update_repulsion(Uold_beta , alpha_coeff, alpha_electrons)

            print '<n_a> = ', occ_alpha, occ_beta

            alpha_hamiltonian = self.__update_hamiltonian(Uold_alpha, Unew_alpha, alpha_hamiltonian)
            beta_hamiltonian = self.__update_hamiltonian(Uold_beta, Unew_beta, beta_hamiltonian)

            print 'RARGH!  WEREWOLVES!'

            Uold_alpha = Unew_alpha
            Uold_beta = Unew_beta

            print '\n'


        print 'Anderson-type Hamiltonians successfully diagonalized.\n'

        self.need_run = False

#        print alpha_coeff
#        print alpha_orbitals
#        print (alpha_coeff[0,0]**2)

#        print alpha_hamiltonian

        return (alpha_orbitals, alpha_coeff, beta_orbitals, beta_coeff)

    def run(self):
        """Calls self.__build_Hamiltonian() to build the alpha and beta spin
           Hamiltonians, and calls self.__scf() to solve the Newns-Anderson
           problem self-consistantly."""
        print '\nNewns.run() issued.'

# Assign number of electrons
        nelectrons = self.keywords['continuum_electrons'] + self.keywords['adsorbate_electrons']
        alpha_electrons, beta_electrons = self.__assign_electrons(nelectrons)
        self.orbital_information['alpha_electrons'] = alpha_electrons
        self.orbital_information['beta_electrons'] = beta_electrons

        print '\nFound %i total electrons:' % nelectrons
        print '%4sAlpha electrons: %i' % (' ', alpha_electrons)
        print '%4sBeta electrons: %i' % (' ', beta_electrons)

        state_range = list(self.keywords['continuum_range'])
        con_ele = int(self.keywords['continuum_electrons'])

        ads_energy = list(self.keywords['adsorbate_energies'])
        ads_ele = int(self.keywords['adsorbate_electrons'])

        vak = float(self.keywords['metal_adsorbate_coupling'])
        Urep = float(self.keywords['U_repulsion'])

        if (state_range) and (con_ele > 2):
            print '\nDiscretizing metal |k> states between %6.3f and %-6.3f eV' % (state_range[0], state_range[1])
            states = self.__assign_electrons(con_ele)
#            print 'Testing:', states
            kstates = self.__discretize(state_range, states[0]-1)
        else:
            print '\nFound 1 metal |k> state with an energy of %6.3f eV' % state_range[0]
            kstates = [state_range[0]]

        self.keywords['kstates'] = kstates

# First build the alpha- and beta-spin Hamiltonians
        print '\nBuilding initial alpha- and beta-spin Anderson-type Hamiltonians.'
        alpha_hamiltonian = self.__build_Hamiltonian(ads_energy, ads_ele,
                                                     kstates, con_ele, vak)
        beta_hamiltonian = self.__build_Hamiltonian(ads_energy, ads_ele,
                                                    kstates, con_ele, vak)

# Right now this will handle a 1-electron adsorbate, will generalize later
# As a 0th-order guess, I assume that <n_{a\sigma}> = 1.  The build_Hamiltonian
# knows nothing about the U term, so we use the self.__update_hamiltonian function
# Since I assume that the spare electron is alpha-spin by convention, we modify
# beta Hamiltonian
        beta_hamiltonian = self.__update_hamiltonian(0.0, Urep, beta_hamiltonian)

        print 'Testing Hamil\'s:'
        print alpha_hamiltonian
        print beta_hamiltonian

        print 'Initial Anderson-type Hamiltonians built.'

# Call the SCF
        alpha_orbitals, alpha_coeff, beta_orbitals, beta_coeff = self.__scf(alpha_hamiltonian, beta_hamiltonian)

# Set the orbital information
        self.orbital_information['alpha_energies'] = alpha_orbitals
        self.orbital_information['alpha_coefficients'] = alpha_coeff

        self.orbital_information['beta_energies'] = beta_orbitals
        self.orbital_information['beta_coefficients'] = beta_coeff

        print 'The following information is now available:\n'
        print '%-30s| %-30s' % ('Information', 'Function')
        print '------------------------------|------------------------------'
        print '%-30s| %-30s' % ('Chemisorption Energy', 'Newns.get_binding_energy()')
        print '\n\n'

    def get_binding_energy(self):
        """Calculates and returns the binding/chemisorption energy"""
        if (self.need_run):
            print 'Hamiltonian not set up or diagonalized.'
            print 'Use Newns.run() to do so.'
            raise RuntimeError('Hamiltonian not setup')

# Aufbau principal:  Sort alpha and beta eigenvalues and then occupy
# In the following, if there is an odd number of electrons, I assume that the
# odd electron is of alpha spin.  In reality, we would have a superposition of 
# spin states, but since I'm treating the orbital energies as degenerate with
# respect to spin, this doesn't matter.
        etot = 0.0 # Total energy, initialized to zero
        e_alpha = sorted(list(self.orbital_information['alpha_energies']))
        alpha_ele = int(self.orbital_information['alpha_electrons'])

        e_beta = sorted(list(self.orbital_information['beta_energies']))
        beta_ele = int(self.orbital_information['beta_electrons'])

# DEBUG:  Checking the memory locations of the new lists
#        print self.orbital_information['alpha_energies']
#        print id(self.orbital_information['alpha_energies'])
#        print e_alpha
#        print id(e_alpha)
#        print '\n'
#        print self.orbital_information['beta_energies']
#        print id(self.orbital_information['beta_energies'])
#        print e_beta
#        print id(e_beta)

# Sum up the total energy, according to spin
        for i in range(alpha_ele):
            etot += e_alpha[i]
        for i in range(beta_ele):
            etot += e_beta[i]

# Get sum of adsorbate states and continuum energies
        ads_energies = list(self.keywords['adsorbate_energies'])
        ads_ele = int(self.keywords['adsorbate_electrons'])
        ads_ele_alpha, ads_ele_beta = self.__assign_electrons(ads_ele)
#        print ads_ele, ads_ele_alpha, ads_ele_beta

        con_energies = list(self.keywords['kstates'])
        con_ele = int(self.keywords['continuum_electrons'])
        con_ele_alpha, con_ele_beta = self.__assign_electrons(con_ele)
#        print con_ele, con_ele_alpha, con_ele_beta

        eads = 0.0
        econ = 0.0

# Here I assume that there isn't any spin splitting occuring, thus each of the
# adsorbate and metal states are degenerate with respect to spin.
        for i in range(ads_ele_alpha):
            eads += ads_energies[i]
        for i in range(ads_ele_beta):
            eads += ads_energies[i]

        for i in range(con_ele_alpha):
            econ += con_energies[i]
        for i in range(con_ele_beta):
            econ += con_energies[i]

# DEBUG: Testing energy values.
#        print 'E_tot =', etot
#        print e_alpha
#        print e_beta
#        print '\nE_ads =', eads
#        print ads_energies
#        print '\nE_con =', econ
#        print con_energies
#        print '\n'

        return etot - (eads + econ)

    def __discretize(self, krange, states):
        """Discretizes a conituum of states, given a range and number of 
        desired states"""
        if (states == 0):
            print '\n0 states does not make sense!'
            raise RuntimeError('Non-sensical states given')

        low = krange[0]
        high = krange[1]
        step = (high - low) / float(states)

        discretize = list()
        for i in range(states + 1):
            discretize.append(low + float(step * i))

        if (i > 1):
            print '%i states were generated.' % states
        else:
            print '%i state was generated.' % states
        return discretize
