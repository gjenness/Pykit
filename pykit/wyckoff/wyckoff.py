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

        self.read_sites(self.keywords['spacegroup'])

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
            print sites

#            print contents[linecount:linecount+int(sites[:-1])/2]
            for line in contents[linecount:linecount+int(sites[:-1])/2]:
                print line.strip()
            linecount += (int(sites[:-1]) / 2 + 1)


#    def read_site_info(self, label, information):
#        print label