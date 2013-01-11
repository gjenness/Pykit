import numpy as np

import os
import warnings

from mastobj import MASTobj


# allowed_keys is a dictionary containing all the allowed keywords for a class
# each value for each key is a tuple, containing an allowed value type, and a 
# default value for that keyword
allowed_keys = {'spacegroup': (int, 1, 'Spacegroup number'),
                'threshold': (float, 0.001, 'Threshold for dectecting sites'),
                'basis': (list, list(), 'List of positions in direct'),
                'print': (bool, False, 'Prints site info to a file'),
               }
            

class Wyckoff(MASTobj):
    def __init__(self, **kwargs):
        MASTobj.__init__(self, allowed_keys, **kwargs)

        self.wyckoff_sites = self.read_sites(self.keywords['spacegroup'])

    def read_sites(self, spacegroup):
        """
        Format of file:

           offset float float float
           sites str
           site name
           site positions
           ....

        Each site position is given on a single line.
        """

        wyckoff_sites = dict()
        coords = ['x', 'y', 'z']

        dirname = os.environ['WYCKOFF_DATA']
        filename = dirname + '/' + str(spacegroup) + '/wyckoff.' + str(spacegroup)
        fsites = open(filename)
        contents = fsites.readlines()
        fsites.close()

        tmp = contents[0].split()
        wyckoff_sites['offset'] = (float(tmp[1]), float(tmp[2]), float(tmp[3]))

        wyckoff_sites['sites'] = contents[1].split()[1:]

        print wyckoff_sites
        linecount = 3
        for sites in wyckoff_sites['sites']:
#            print sites
            siteslist = list()
#            print contents[linecount:linecount+int(sites[:-1])/2]
            for line in contents[linecount:linecount+int(sites[:-1])/2]:
                tmp = line.strip().split()
                tmp1 = list()
                for tm in tmp:
                    if ('+' in tm):
                       tm = tm.split('+')
                       if ('x' in tm[0]) or ('y' in tm[0]) or ('z' in tm[0]):
                           tm[1] = float(tm[1])
                           tm = (tm[0], tm[1])
                    else:
                        try:
                            tm = float(tm)
                        except ValueError:
                            pass
                    tmp1.append(tm)
#                print tmp1
                siteslist.append(tmp1)

            wyckoff_sites[sites] = siteslist
            linecount += (int(sites[:-1]) / 2 + 1)
#        print wyckoff_sites

        return wyckoff_sites

    def get_wyckoff_sites(self, positions):
        print self.wyckoff_sites['sites']
        sites = list(self.wyckoff_sites['sites'])
        site_index = dict()

        check = list(self.wyckoff_sites[sites[0]])
        offset = np.array(self.wyckoff_sites['offset'])
        print check
        print offset

# Create a list of the lower symmetry sites --- these depend only on fixed
# coordinates, so it's best to check for these first!
# Expand out with the offset, if a value is equal to 1, set it to 0 due to PBCs
        for che in self.wyckoff_sites[sites[0]]:
           new = np.array(che) + offset
           for i in range(3):
               if new[i] == 1.0:
                   new[i] = 0.0
           check.append(new)

        check = np.array(check)
        print self.keywords['threshold']
        for i in range(len(positions)):
            tmp = list()
            count = 1
#            print 'Atom %i, %10.2f%10.2f%10.2f' % (i, abs(positions[i, 0] - che[0]),
#                                          abs(positions[i, 1] - che[1]),
#                                          abs(positions[i, 2] - che[2]) )
            for che in check:
#                print i, positions[i], che, (positions[i] - che)
                if (abs(positions[i, 0] - che[0]) < self.keywords['threshold']) and \
                   (abs(positions[i, 1] - che[1]) < self.keywords['threshold']) and \
                   (abs(positions[i, 2] - che[2]) < self.keywords['threshold']):
                    print 'Found a %s site! Atom:  %i' % (sites[0], i)
                    count += 1
