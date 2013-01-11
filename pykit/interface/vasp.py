from pykit.interface.program import Program

vasp_keywords = {'potim': (int, 'default', 'MD timestep'),
                 'encut': (float, 'default', 'Planewave energy cutoff (eV)'),
                 'ibrion': (int, 2, 'Optimization method, default conjugant-gradient')
                }

class Vasp(Program):
    def __init__(self, **kwargs):
        Program.__init__(self, **kwargs)

        self.vasp_keys = self.derived_keys.copy()
