from UserDict import IterUserDict # May decide to not need this!

from pykit.util.mastobj import MASTobj

program_keywords = {'subdirectory': (str, None, 'Creates a subdirectory to run a job in'),
                   }

class Program(MASTobj):
    """Generic interface to a DFT program."""
    def __init__(self, **kwargs):
        MASTobj.__init__(self, **kwargs)

# self.program_data will hold any data that is retrieved from the output
        self.program_data = dict()
        self.program_keys, self.derived_keys = self.split_keywords(**kwargs)

# General utility functions related to running the program
    def split_keywords(self, **kwargs):
        """Takes kwargs and breaks them into general keywords related to
           the interface, and those related to the program. Returns them as a
           tuple"""
        progamkeys = dict()
        derivedkeys = dict()

        for key, value in kwargs.items():
            if key in program_keywords:
                programkeys[key] = value
            else:
                derivedkeys[key] = value

         return (programkeys, derivedkeys)

# Functions the are uniquely defined for each program.
    def _write_input(self):
        """Define in derived class"""
        raise NotImplementedError

    def _read_output(self):
        """Define in derived class"""
        raise NotImplementedError

    def _run(self):
        """Define in derived class"""
        raise NotImplementedError

#  Anything below this should not need to be redefined!
    def get_energy(self):
        """Returns the energy"""
        raise NotImplementedError

    def get_unit_cell(self):
        """Returns the unit cell"""
        raise NotImplementedError

    def get_coordinates(self):
        """Returns the coordinates"""
        raise NotImplementedError

    def get_forces(self):
        """Returns the forces"""
        raise NotImplementedError

    def get_stresses(self)
        """Returns the stresses"""
        raise NotImplementedError


