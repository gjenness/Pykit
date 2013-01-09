import numpy as np
import numpy.linalg as la


class Newns:
    def __init__(self, adsorbate=list(), srange=tuple(), prange=tuple(), drange=tuple(), states=0):
        print 'Newns-Anderson solver'

        self.adsorbate = adsorbate
        self.srange = srange
        self.prange = prange
        self.drange = drange
        self.states = states

        self.need_run = False

    def __build_Hamiltonian(self):
        """Builds the Newns-Anderson Hamiltonian"""

        print '\nBuilding Newns-Anderson Hamiltonian'
        print 'Discretizing metal |k> states'
        if (self.srange):
            self.sstates = self.__discretize(self.srange, self.states)
        elif (not self.srange):
            self.sstates = list()

        if (self.prange):
            self.pstates = self.__discretize(self.prange, self.states)
        elif (not self.prange):
            self.pstates = list()

        if (self.drange):
            self.dstates = self.__discretize(self.drange, self.states)
        elif (not self.drange):
            self.dstates = list()

        # Length of various lists
        las = len(self.adsorbate)
        lss = len(self.sstates)
        lps = len(self.pstates)
        lds = len(self.dstates)

        nelements = las + lss + lps + lds
        hamiltonian = np.zeros( (nelements, nelements) )

# For debugging purposes
#        print len(sstates)
#        print sstates
#        print len(pstates)
#        print len(dstates)

        # Rock that diagonal!
        # Some counters
        nas = 0
        nss = 0
        nps = 0
        nds = 0


        for i in range(nelements):
            if (i < las):
                hamiltonian[i, i] = self.adsorbate[nas]
                nas += 1
            elif (i < (las + lss)):
                hamiltonian[i, i] = self.sstates[nss]
                nss += 1
            elif (i < (las + lss + lps)):
                hamiltonian[i, i] = self.pstates[nps]
                nps += 1
            elif (i < (las + lss + lps + lds)):
                hamiltonian[i, i] = self.dstates[nds]
                nds += 1

        # Rock that off-diagonal!
        for i in range(las, nelements):
            for j in range(las):
                if (i != j):
                    hamiltonian[i, j] = -9.5e0

                hamiltonian[j, i] = hamiltonian[i, j]

        print 'Hamiltonian built!\n'
# For debugging purposes
        txt = open('hamil.txt', 'w')
        for i in range(nelements):
            for j in range(nelements):
                txt.write('%9.3f' % hamiltonian[i, j])
            txt.write('\n')
        txt.close()

        return hamiltonian

    def run(self):
        """Calls self.build_Hamiltonian, and diagonalizes it"""
        # First build the Hamiltonian
        hamiltonian = self.__build_Hamiltonian()
        print 'Diagonalizing Hamiltonian\n'

        self.eigenvalues, self.eigenvectors = la.eig(hamiltonian)
        print 'Hamiltonian successfully diagonalized\n'
#        print self.eigenvalues

        self.need_run = True

    def get_binding_energy(self):
        if (not self.need_run):
            self.run()

        etot = sum(self.eigenvalues)
        eads = sum(self.adsorbate)
        eks  = sum(self.sstates) + sum(self.pstates) + sum(self.dstates)

#        print self.eigenvalues
        print 'Separated', etot, eads, eks
        print '2 times:', etot - 2 * (eads + eks)
        print '1 times:', etot - (eads + eks)

        return etot - 2 * (eads + eks)

    def __discretize(self, krange, states):
        """Discretizes a conituum of states, given a range and number of 
        desired states"""
        low = krange[0]
        high = krange[1]
        step = (high - low) / float(states)

        discretize = list()
        for i in range(states + 1):
            discretize.append(low + float(step * i))

        return discretize
