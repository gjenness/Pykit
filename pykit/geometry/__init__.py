import os
import numpy as np


def geomread(filename, format=None):
    if ('poscar' in filename.lower()):
        from pykit.geometry.vasp import Poscar
        return Poscar(filename)
    elif (filename.lower().endswith('.xyz')):
        from pykit.geometry.tinker import Tinker
        return Tinker(filename)
    else:
        print 'File type', filename, 'not recognized.'
        raise RuntimeError('Bad file type.')

def geomconvert(files, formats=[]):
    file1 = files[0]
    file2 = files[1]

    geom1 = geomread(file1)
    unit_cell = geom1.get_unit_cell()
    positions = geom1.get_positions()
    atom_list = geom1.get_atom_list()

    if ('poscar' in file2.lower()):
        from pykit.geometry.vasp import Poscar
        geom2 = Poscar(file2, status=None)
    elif ('xyz' in file2.lower()):
        from pykit.geometry.tinker import Tinker
        geom2 = Tinker(file2, status=None)
    else:
        print 'File type', filename, 'not recognized.'
        raise RuntimeError('Bad file type.')

    geom2.set_unit_cell(unit_cell)
    geom2.set_positions(positions)
    geom2.set_atom_list(atom_list)

    geom2.write_posfile(file2)
 

