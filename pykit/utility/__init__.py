"""
Utility functions for PyKit
"""


def convert2list(string):
    """Checks a string to see if it's a list or not.
    If not a list, converts it to a list"""
    if isinstance(string, list):
        return string
    else:
        return [string]

def integrate(grid=None, data=None):
    """Performs a numerical integration using a non-uniform trapezoid
    approximation"""
    if (len(grid) != len(data)):
        print 'Grid length: ', len(grid)
        print 'Data length: ', len(data)
        raise RuntimeError('Mismatch in data length')

    ngrid = len(grid)
#    print 'Found', ngrid, 'grid points.'
#    print 'Integration end-points are: ', grid[0], grid[-1]

    kernal = 0.0
    for i in range(ngrid-1):
#        print '\n'
#        print grid[i+1], grid[i]
#        print data[i+1], data[i]
        kernal += (grid[i+1] - grid[i]) * (data[i+1] + data[i])

    return kernal * 0.5
