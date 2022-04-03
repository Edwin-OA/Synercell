import sys
import os
from reaktoro import *
import math

this_dir = os.path.dirname(__file__)


class BioMolecule:
    """

    Simplified version of BioMolecule from NutMEG's synthesis submodule

    """
    name = ''
    Mr = 0. # molecular mass
    conc_mol_per_l = None	# mol/l
    conc_mol_per_cell = None # mol/cell

    std_formation_gibbs = None # kJ/mol
    std_formation_R = None #kJ/mol

    stdbio_formation_gibbs = None # (biological standard) kJ/mol
    stdbio_formation_R = None # (biological standard) kJ/mol


    def __init__(self, name, Mr, N, T=298.15, P=101325.,
      pH=7., z=0., I=0. ):
        self.name = name
        self.Mr = Mr

        self.GetThermoParams(T, P, N, pH=pH, I=I, z=z)


    def GetThermoParams(self, T, P, N, pH=7., I=0., z=0.):

        db = Database(this_dir+"/supcrt07-organics_AABB.xml")
        thermodynamics = Thermo(db)
        self.std_formation_gibbs = thermodynamics.standardPartialMolarGibbsEnergy(T, P, self.name).val *0.001

        AABB = thermodynamics.standardPartialMolarGibbsEnergy(T, P, 'AABB').val *0.001

        self.std_formation_R = self.std_formation_gibbs - AABB

        self.stdbio_formation_gibbs = BioMolecule.biogibbs(
          self.std_formation_gibbs, T, pH, N, I=I, z=z)

        stdbio_AABB = BioMolecule.biogibbs(
          AABB, T, pH, 4, I=I, z=z)

        self.stdbio_formation_R = BioMolecule.biogibbs(
          self.std_formation_R, T, pH, N-4, I=I, z=0.)

    @staticmethod
    def biogibbs(stdGibbs, T, pH, N, I=0., z=0):
        """Uses eqs from Alberty 1998 'Calculation  of  Standard Transformed
        Gibbs Energiesand Standard Transformed Enthalpies of Biochemical
        Reactants', Archives of Biochemistry and Biophysics
        """
        if I == 0.:
            return stdGibbs - (N*0.00831*T*math.log(10**(-pH)))
        else:
            return stdGibbs - (2.91482*z*z*(I**0.5)/(1+(B*(I**0.5)))) - N*(
              - 2.91482*(I**0.5)/(1+(B*(I**0.5)))
              + 0.00831*T*math.log(10**(-pH)))
