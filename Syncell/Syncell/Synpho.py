"""
Synpho, a class for calculating the energetic costs assoicated with
synthesising a cell membrane.

An instance describes the membrane of a cell as a set and of its
individual parts at selected environmental conditions. At present, the
only condition which varies behaviour is temperature.

"""

import sys,os, math, statistics, re
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.Seq import transcribe
from BioMolecule import BioMolecule

from GCAtools import *
from collections import Counter
from copy import deepcopy

class Synpho:
    """
    Class for calculating the energetic costs assoicated with
    synthesising a membrane.

    Properties
    ----------
    T : float
        temperature governing themodynamic constants.
    Metabolites : dict
        dictionary of thermodynamic data for metabolites.
    cell : Syncell.Cell
        instance of Cell, collecting microbial properties.
    """

    def __init__(self, T, cell):
        self.T = T
        self.Metabolites = Synpho.get_Metabolites_dict(T)
        self.cell = cell

    def update_T(self, T):
        """
        Update the properties of this Synpho object relating to new
        temperature T.
        """
        self.T = T
        self.Metabolites = Synpho.get_Metabolites_dict(T)

    @staticmethod
    def get_Metabolites_dict (T):
        """
        Return a ditionary of the relevant metabolites, along with
        relevant 'backbone' species.

        The dictionary contains lists with the following indexes:
            0= Molecular formula
            1= Molecular weight (Da)
            2= dGf
            3= E. coli - absolute intracellular concentration (mol/L)
                (Bennett et al. 2009 & Park et al. 2016)
            4= Yeast - absolute intracellular concentration (mol/L)
                (Bennett et al. 2009 & Park et al. 2016)
            5= Mammalian - absolute intracellular concentration (mol/L)
                (Bennett et al. 2009 & Park et al. 2016)

        """
        _dAS = BioMolecule('Deoxyadenosine(aq)', 251.24192, 13, T=T)
        _A = BioMolecule('Adenosine(aq)', 267.24132, 13, T=T)

        # sugars
        _rib = BioMolecule('Ribose(aq)', 150.1299, 10, T=T)
        _drib = BioMolecule('Deoxyribose(aq)', 134.1305, 10, T=T)
        _R5P = BioMolecule('Ribose-5-Phosphate-2', 228.093921, 9, T=T, z=-2)

        # Phosphates
        _H2PO4 = BioMolecule('H2PO4-', 96.987241, 2, T=T, z=-1)
        _H3PO4 = BioMolecule('H3PO4(aq)', 97.995181, 3, T=T)
        _HPO4 = BioMolecule('HPO4--', 95.979301, 1, T=T, z=-2)
        _PO4 = BioMolecule('PO4---', 94.971361, 0, T=T, z=-3)
        _P2O7 = BioMolecule('P2O7----', 173.943322, 0, T=T, z=-4)
        # _ribphos = BioMolecule('Ribose-5-Phosphate-2', 150.1299, 10, T=T, z=-2)

        #Amino acids
        _C5H10N2O3 = BioMolecule('Glutamine(aq)', 146.1445, 10, T=T)
        _C3H7NO3 = BioMolecule('Serine(aq)', 105.09258, 7, T=T)

        # Other Metabolites
        _Water = BioMolecule('H2O(l)', 18, 2, T=T)
        _H2 = BioMolecule('H2(g)', 2.01, 2, T=T)
        _NH4 = BioMolecule('NH4+', 18.03, 4, T=T, z=1)
        _C6H12O6 = BioMolecule('Glucose(aq)', 180.15588, 12, T=T)
        _C3H2O4 = BioMolecule('Malonate--', 102.04558, 2, T=T, z=-2)
        _C10H16N5O13P3 = BioMolecule('+H4ATP-(aq)', 507.181023, 16, T=T)
        _C10H15N5O10P2 = BioMolecule('+H3ADP1-(aq)', 427.201122, 15, T=T)
        _O2 = BioMolecule('O2(aq)', 31.9988, 0, T=T)
        _CO2 = BioMolecule('CO2(aq)', 44.0095, 0, T=T)

        # Bonds (KJ/mol)
        Phosphate = _HPO4.stdbio_formation_gibbs
        dGf_PhdB = 22.175 # (Reversed) Phosphodiester bond At 25 °C, pH 7  Dickson K. 2000
        ester_bond = _P2O7.stdbio_formation_gibbs - ((2*_HPO4.stdbio_formation_gibbs) - _Water.stdbio_formation_gibbs)
        #ester_bond = (_rib.stdbio_formation_gibbs + Phosphate) - (_R5P.stdbio_formation_gibbs - _Water.stdbio_formation_gibbs)
        N_glycosidic_bond = _dAS.stdbio_formation_gibbs - ((_A.stdbio_formation_gibbs + _drib.stdbio_formation_gibbs) - _Water.stdbio_formation_gibbs)
        _OH = _rib.stdbio_formation_gibbs - _drib.stdbio_formation_gibbs

        #Backbones
        #DNA phosphate backbone = Phosphate + ester bond + Deoxyribose
        DNA_phosphate_BB = (
        Phosphate + ester_bond + _drib.stdbio_formation_gibbs - _Water.stdbio_formation_gibbs)
        #RNA phosphate backbone = Phosphate + ester bond + Ribose
        RNA_phosphate_BB=  (
        Phosphate + ester_bond + _rib.stdbio_formation_gibbs - _Water.stdbio_formation_gibbs)
        #Phosohate + ester bond
        phosphoester = (
        (Phosphate + ester_bond) - _Water.stdbio_formation_gibbs)

        #Phospholipids metabolites
        #Metabolites energies
        dGf_H = 0
        w_H = 1
        dGf_e = 0
        '''dGf'''
        glucose = _C6H12O6.stdbio_formation_gibbs
        ATP = _C10H16N5O13P3.stdbio_formation_gibbs
        ADP = _C10H15N5O10P2.stdbio_formation_gibbs
        glutamine = _C5H10N2O3.stdbio_formation_gibbs
        serine = _C3H7NO3.stdbio_formation_gibbs
        O2 = _O2.stdbio_formation_gibbs
        CO2 = _CO2.stdbio_formation_gibbs
        H2O = _Water.stdbio_formation_gibbs
        PO4 = _PO4.stdbio_formation_gibbs
        H2PO4 = _H2PO4.stdbio_formation_gibbs
        Phosphate = _HPO4.stdbio_formation_gibbs
        malonate = _C3H2O4.stdbio_formation_gibbs
        dGf_R5P = _R5P.stdbio_formation_gibbs
        '''Weight'''
        w_glucose = _C6H12O6.Mr
        w_ATP = _C10H16N5O13P3.Mr
        w_ADP = _C10H15N5O10P2.Mr
        w_glutamine = _C5H10N2O3.Mr
        w_serine = _C3H7NO3.Mr
        w_O2 = _O2.Mr
        w_CO2 = _CO2.Mr
        w_H2O = _Water.Mr
        w_PO4 = _PO4.Mr
        w_H2PO4 = _H2PO4.Mr
        w_Phosphate = _HPO4.Mr
        w_malonate = _C3H2O4.Mr
        w_dGf_R5P = _R5P.Mr

        '''dGf by group contribution'''
        glucose_6_p = ((ATP - ADP) + glucose)
        G3P = (((glucose * 0.5) + (ATP * 1)) - ((ADP * 1)))
        choline = (((serine * 3) + (ATP * 7.4) + (H2O * 0.2)) - ((ADP * 7.8) + (PO4 * 6.6) + (dGf_H * 8.8)))
        glycerol = (((glucose_6_p * 0.5) + (ATP * 2.5) + (O2 * 1.5)) - ((ADP * 2.5) + (PO4 * 3)))
        pyruvate = (((G3P * 1) + (ADP * 2 ) + (PO4*1)) - ((ATP * 2) + (H2O * 1)))
        palmitate = (((pyruvate * 8) + (malonate*2)) - ((CO2 * 14) + (H2O * 2)))
        oleate = (((pyruvate * 9) + (malonate*2)) - ((CO2 * 15) + (H2O * 3)))
        POPC = (oleate + palmitate + choline + glycerol + H2PO4) - (4 * H2O)
        # POPC = (((glucose * 5) + (serine * 1) + (ATP * 1) + (pyruvate*10) + (malonate * 1)) - ((ADP * 1) + (H2O * 14) + (CO2 * 24)))

        '''weight by group contribution'''
        w_glucose_6_p = ((w_ATP - w_ADP) + w_glucose)
        w_G3P = (((w_glucose * 0.5) + (w_ATP * 1)) - ((w_ADP * 1)))
        w_choline = (((w_serine * 3) + (w_ATP * 7.4) + (w_H2O * 0.2)) - ((w_ADP * 7.8) + (w_PO4 * 6.6) + (w_H * 8.8)))
        w_glycerol = (((w_glucose_6_p * 0.5) + (w_ATP * 2.5) + (w_O2 * 1.5)) - ((w_ADP * 2.5) + (w_PO4 * 3)))
        w_pyruvate = (((w_G3P * 1) + (w_ADP * 2 ) + (w_PO4*1)) - ((w_ATP * 2) + (w_H2O * 1)))
        w_palmitate = (((w_pyruvate * 8) + (w_malonate*2)) - ((w_CO2 * 14) + (w_H2O * 2)))
        w_oleate = (((w_pyruvate * 9) + (w_malonate*2)) - ((w_CO2 * 15) + (w_H2O * 3)))
        #Phospholipids
        # w_POPC = (((w_glucose * 5) + (w_serine * 1) + (w_ATP * 1) + (w_pyruvate * 10) + (w_malonate * 1)) - ((w_ADP * 1) + (w_H2O * 14) + (w_CO2 * 24)))
        w_POPC = (w_oleate + w_palmitate + w_choline + w_glycerol + w_H2PO4) - (4 * w_H2O)


        Metabolites = {
        'glucose': {'name':'glucose', 'MF': 'C6H12O6', 'type': 'metabolite', 'MW': _C6H12O6.Mr, 'dGf': _C6H12O6.stdbio_formation_gibbs, 'EC[]': 2.50E-05,'SC[]': 5.31E-03,'M[]': 6.75E-04},
        'ATP': {'name': 'ATP', 'MF': 'C10H16N5O13P3', 'type': 'metabolite', 'MW': _C10H16N5O13P3.Mr, 'dGf': _C10H16N5O13P3.stdbio_formation_gibbs, 'EC[]': 9.63E-03,'SC[]': 1.93E-03,'M[]': 4.67E-03},
        'ADP': {'name': 'ADP', 'MF': 'C10H15N5O10P2', 'type': 'metabolite', 'MW': _C10H15N5O10P2.Mr, 'dGf': _C10H15N5O10P2.stdbio_formation_gibbs, 'EC[]': 5.55E-04,'SC[]': 4.88E-04,'M[]': 5.69E-04},
        'glutamite': {'name': 'glutamite', 'MF': 'C5H10N2O3', 'type': 'amino acid', 'MW': _C5H10N2O3.Mr, 'dGf': _C5H10N2O3.stdbio_formation_gibbs, 'EC[]': 3.81E-03,'SC[]': 3.55E-02,'M[]': 1.62E-02},
        'serine': {'name': 'serine', 'MF': 'C3H7NO3', 'type': 'amino acid', 'MW': _C3H7NO3.Mr, 'dGf': _C3H7NO3.stdbio_formation_gibbs, 'EC[]': 6.8E-5,'SC[]': 6.8E-5,'M[]': 1.81E-5},
        'O2': {'name': 'O2', 'MF': 'O2', 'type': 'metabolite', 'MW': _O2.Mr, 'dGf': _O2.stdbio_formation_gibbs, 'EC[]': '?','SC[]':'?','M[]':'?'},
        'CO2': {'name': 'CO2', 'MF': 'CO2', 'type': 'metabolite', 'MW': _CO2.Mr, 'dGf': _CO2.stdbio_formation_gibbs, 'EC[]': 7.52E-05,'SC[]': 8.16E-05,'M[]': 7.63E-03},
        'H2O': {'name': 'H2O', 'MF': 'H2O', 'type': 'metabolite', 'MW': _Water.Mr, 'dGf': _Water.stdbio_formation_gibbs, 'EC[]': 5.56E-14,'SC[]': 6.11E-12,'M[]': 1.668E-10},
        'PO4': {'name': 'PO4', 'MF': 'PO4', 'type': 'metabolite', 'MW': _PO4.Mr, 'dGf': _PO4.stdbio_formation_gibbs, 'EC[]': 2.39E-02,'SC[]': 4.93E-02,'M[]': 5.83E-03},
        'H2PO4': {'name': 'H2PO4', 'MF': 'H2PO4', 'type': 'metabolite', 'MW': _H2PO4.Mr, 'dGf': _H2PO4.stdbio_formation_gibbs, 'EC[]': 2.39E-02,'SC[]': 4.93E-02,'M[]': 5.83E-03},
        'phosphate': {'name':'phosphate', 'MF': 'HPO4', 'type': 'metabolite', 'MW': _HPO4.Mr, 'dGf': _HPO4.stdbio_formation_gibbs, 'EC[]': 2.39E-02,'SC[]': 4.93E-02,'M[]': 5.83E-03},
        'malonate': {'name': 'malonate', 'MF': 'C3H2O4', 'type': 'metabolite', 'MW': _C3H2O4.Mr, 'dGf': _C3H2O4.stdbio_formation_gibbs, 'EC[]': 7.26E-05,'SC[]': 7.26E-05,'M[]': 7.26E-05},
        'choline': {'name': 'choline', 'MF': 'C5H14NO', 'type': 'metabolite', 'MW': w_choline, 'dGf': choline, 'EC[]': 3.28E-11,'SC[]': 3.28E-114,'M[]': 3.28E-11},
        'glycerol': {'name': 'glycerol', 'MF': 'C3H8O3', 'type': 'metabolite', 'MW': w_glycerol, 'dGf': glycerol, 'EC[]': 1.61E-04,'SC[]': 2.81E-04,'M[]': 4.90E-05},
        'pyruvate': {'name': 'pyruvate', 'MF': 'C3H4O3', 'type': 'metabolite', 'MW': w_pyruvate, 'dGf': pyruvate, 'EC[]': 3.66E-03,'SC[]': 9.40E-03,'M[]': 5.88E-03},
        'palmitate': {'name': 'palmitate', 'MF': 'C16H32O2', 'type': 'metabolite', 'MW': w_palmitate, 'dGf': palmitate, 'EC[]': 1.12E-04,'SC[]': 1.12E-04,'M[]': 1.12E-04},
        'oleate': {'name': 'oleate', 'MF': 'C18H34O2', 'type': 'metabolite', 'MW': w_oleate, 'dGf': oleate, 'EC[]': 6.22E-05,'SC[]': 6.22E-05,'M[]': 6.22E-055},
        'POPC': {'name': 'POPC', 'MF': 'C42H82NO8P', 'type': 'phospholipid', 'MW': w_POPC , 'dGf': POPC, 'EC[]': '?','SC[]': '?','M[]': '?'},
        'deoxyribose': {'name': 'deoxyribose', 'MF': 'C5H10O4', 'type': 'metabolite', 'MW': _drib.Mr, 'dGf': _drib.stdbio_formation_gibbs, 'EC[]': 3.03E-4,'SC[]': '?','M[]': '?'},
        'ribose': {'name': 'ribose', 'MF': 'C5H10O5', 'type': 'metabolite', 'MW': _rib.Mr, 'dGf': _rib.stdbio_formation_gibbs, 'EC[]': 7.83E-5,'SC[]': 1.52E-4,'M[]': '?'},
        'ester_bond': {'name': 'ester_bond', 'MF': 'RCO2R', 'type': 'bond', 'MW': 'none', 'dGf': ester_bond, 'EC[]': 'none','SC[]': 'none','M[]': 'none'},
        'ROH': {'name': 'glycosidic_bond', 'MF': 'ROH', 'type': 'bond', 'MW': 'none', 'dGf': N_glycosidic_bond, 'EC[]': 'none','SC[]': 'none','M[]': 'none'},
        'Phosphate + ester bond': {'name':'phosphoester', 'MF': 'Phosphate + ester bond', 'type': 'bond', 'MW': 'none', 'dGf': phosphoester, 'EC[]': 'none','SC[]': 'none','M[]': 'none'},
        'Phosphate + ester bond + Deoxyribose': {'name':'DNA_PBB', 'MF': 'Phosphate + ester bond + Deoxyribose', 'type': 'backbone', 'MW': 'none', 'dGf': DNA_phosphate_BB, 'EC[]': 'none','SC[]': 'none','M[]': 'none'},
        'Phosphate + ester bond + Ribose': {'name': 'RNA_PBB', 'MF': 'Phosphate + ester bond + Ribose', 'type': 'backbone', 'MW': 'none', 'dGf': RNA_phosphate_BB, 'EC[]': 'none','SC[]': 'none','M[]': 'none'},
        'OH': {'name': 'OH', 'MF': 'OH', 'type': 'functional group', 'MW': 'none', 'dGf': _OH, 'EC[]': 'none','SC[]': 'none','M[]': 'none'}}


        return Metabolites

    '''Weight of the biomolecule'''
    @staticmethod
    def get_weight_p(dictionary):
        tot_weight = 0
        for pl in dictionary:
            if dictionary[pl]['count'] != 0:
                P_weight = dictionary[pl]['MW'] * dictionary [pl]['count']
                tot_weight += P_weight
        return tot_weight

    '''Energy of the membrane'''
    ''' Formation energy'''
    @staticmethod
    def get_dGf_membrane (dictionary):
        membrane_dGf = 0
        for pl in dictionary:
            if dictionary[pl]['count'] != 0:
                pl_dGf = dictionary[pl]['dGf'] * dictionary [pl]['count']
                membrane_dGf += pl_dGf

        return membrane_dGf

    ''' Reaction energy'''
    @staticmethod
    def get_dGr(dGf_product, dGf_reactant):
        dGr = dGf_product - dGf_reactant
        return dGr

    '''Reaction quotient '''
    def get_lnQ_pl(self, PC, N_biomolecules, dictionary):
        lnQ = math.log(PC / N_biomolecules)
        for pl in dictionary:
            if dictionary[pl]['count'] != 0:
                if self.cell.get_model() == "EC":
                    lnQ -= math.log(dictionary['glucose']['EC[]']) + math.log(dictionary['serine']['EC[]']) + math.log(dictionary['ATP']['EC[]']) + math.log(dictionary['pyruvate']['EC[]']) + math.log(dictionary['malonate']['EC[]'])
                elif self.cell.get_model() == "SC":
                    lnQ -= math.log(dictionary['glucose']['SC[]']) + math.log(dictionary['serine']['SC[]']) + math.log(dictionary['ATP']['SC[]']) + math.log(dictionary['pyruvate']['SC[]']) + math.log(dictionary['malonate']['SC[]'])
                elif self.cell.get_model() == "M":
                    lnQ -= math.log(dictionary['glucose']['M[]']) + math.log(dictionary['serine']['M[]']) + math.log(dictionary['ATP']['M[]']) + math.log(dictionary['pyruvate']['M[]']) + math.log(dictionary['malonate']['M[]'])

        return lnQ


    def membrane_synthesis(self):
        """
        Calculate and return the free energy of synthesising one dry g of the
        membrane of this instance's Cell.
        """

        '''Reading block'''
        phospholipid_list = []
        for m in self.Metabolites:
            if self.Metabolites[m]['type'] == 'phospholipid':
                phospholipid_list.append(self.Metabolites[m]['name'])
        '''Energy of the membrane'''
        pl_dGf_list = []
        pl_dGr_list = []
        pl_reactionGs = []
        for pl in phospholipid_list:
            # Phospholipid diversity count and length
            pl_cur_dict = deepcopy(self.Metabolites)
            if pl in pl_cur_dict[pl]['name']:
                pl_n_count = Counter(phospholipid_list)

        for pl, pl_dictionary in pl_cur_dict.items():
            pl_dictionary['count'] = pl_n_count.get(pl ,0) #Add 1 to the dictionary counter if the PL is present

        '''Weight of the phospholipid sample'''
        membrane_weight = Synpho.get_weight_p(pl_cur_dict)
        # print('Weight of the phospholipid', membrane_weight)

        '''Membrane dGf '''
        dGf_membrane = Synpho.get_dGf_membrane (pl_cur_dict)
        pl_dGf_list.append(dGf_membrane)
        # print('Formation energy of the membrane', dGf_membrane, 'KJ/mol at', T, 'Kelvin')

        '''Membrane dGr '''
        pl_product = dGf_membrane + (self.Metabolites['ADP']['dGf'] + (self.Metabolites['H2O']['dGf'] * 14) + (self.Metabolites['CO2']['dGf'] * 24))
        pl_reactant = (self.Metabolites['glucose']['dGf'] * 5) + self.Metabolites['serine']['dGf'] + self.Metabolites['ATP']['dGf'] + (self.Metabolites['pyruvate']['dGf'] * 10) + self.Metabolites['malonate']['dGf']
        dGr_membrane = Synpho.get_dGr(pl_product, pl_reactant)
        pl_dGr_list.append(dGr_membrane)
        # print('Reaction energy of the membrane', dGr_membrane, 'KJ/mol at', T, 'Kelvin')

        '''Molarity'''
        pl_N = (self.cell.get_lipids_weight() / self.Metabolites['POPC']['MW']) * 6.0221e23 # Avogadro constant
        av_pl_w = membrane_weight / len(phospholipid_list)
        POPC_moles = (self.cell.get_lipids_weight() / self.Metabolites['POPC']['MW'])
        membrane_molarity =  POPC_moles / self.cell.get_volume()
        # print('membrane_molarity',membrane_molarity)

        '''Reaction quotient'''
        R = 0.008314472 #Gas constant KJ mol K
        for pl, pl_dGr in zip(phospholipid_list, pl_dGr_list):
            pl_lnQ = self.get_lnQ_pl(membrane_molarity, pl_N, pl_cur_dict)
            # print('pl_lnQ', pl_lnQ)

            '''Molar Gibbs Energy'''
            pl_dG = (pl_dGr + (R * self.T * pl_lnQ))
            pl_reactionGs.append(pl_dG)
        # print('Molar gibbs energy', pl_dG)

        '''Membrane energy (KJ/mol)'''
        membrane_energy = get_omic_energy(pl_reactionGs, av_pl_w)

        return membrane_energy
