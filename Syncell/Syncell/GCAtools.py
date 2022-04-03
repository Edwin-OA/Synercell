'''Last edit 02.03.2022 - Pete'''
"""
GCAtools.py

This file contains useful static methods which are useful throughout
Syncell's classes.
"""


from BioMolecule import BioMolecule
from statistics import mean
import math, statistics

'''Weight of the biomolecule'''
def get_weight(dictionary):
    biomolecule_w = 0
    for BB in dictionary:
        BB_weight = dictionary[BB]['MW'] * dictionary [BB]['count']
        biomolecule_w += BB_weight
    return biomolecule_w


'''Reaction Quotient'''
def get_lnQ(PC, biomolecule, N_biomolecules, dictionary, celltype):
    lnQ = math.log(PC / N_biomolecules)
    for BB in biomolecule:
        if BB in dictionary.keys():
                                if celltype == "EC":
                                    lnQ -= math.log(dictionary[BB]['EC[]'])
                                elif celltype == "SC":
                                    lnQ -= math.log(dictionary[BB]['SC[]'])
                                elif celltype == "M":
                                    lnQ -= math.log(dictionary[BB]['M[]'])
    return lnQ

'''Alternative Reaction Quotient'''
def get_alt_lnQ(PC, biomolecule, N_biomolecules, dictionary, celltype):
    bonds = len(biomolecule)

    # this dictionary contains the stoichiomety of the formation reactions
    # of the nucleotides
    NT_contribs = {
      'A': {
        'Glucose':0.5, 'ATP':1.,'Glutamine': 2.5,'O2':6.5, 'nuc':1., 'ADP':1.,'Water':9.5,'CO2':5.5},
      'C' : {
        'Glucose':0.5, 'ATP':1.,'Glutamine': 1.5,'O2':2,'nuc':1., 'ADP':1.,'Water':4.5,'CO2':1.5},
      'G': {
        'Glucose':0.5, 'ATP':1.,'Glutamine': 2.5,'O2':7.,'nuc':1., 'ADP':1.,'Water':9.5,'CO2':5.5},
      'T': {
        'Glucose':0.5, 'ATP':3.,'Glutamine': 3.5,'O2':2.,'nuc':1., 'ADP':4.,'Water':7.5,'CO2':0.5}
      }

    # begin the log-quotient with the DNA molecule
    lnQ = math.log(PC/N_biomolecules)

    for BB in biomolecule:
        if BB in dictionary.keys():

            # add the rest of the products so they are on the top of Q
            if celltype == "EC":
                lnQ += math.log((dictionary['ADP']['EC[]'] * NT_cotribs[BB]['ADP']) + (dictionary['CO2']['EC[]'] * NT_cotribs[BB]['CO2']) + (dictionary['H2O']['EC[]'] * NT_cotribs[BB]['Water']))
            if celltype == "SC":
                lnQ += math.log((dictionary['ADP']['SC[]'] * NT_cotribs[BB]['ADP']) + (dictionary['CO2']['SC[]'] * NT_cotribs[BB]['CO2']) + (dictionary['H2O']['SC[]'] * NT_cotribs[BB]['Water']))
            if celltype == "M":
                lnQ += math.log((dictionary['ADP']['M[]'] * NT_cotribs[BB]['ADP']) + (dictionary['CO2']['M[]'] * NT_cotribs[BB]['CO2']) + (dictionary['H2O']['M[]'] * NT_cotribs[BB]['Water']))

            # take away the reactants so they are on the bottom of Q
            if celltype == "EC":
                lnQ -= math.log((dictionary['ATP']['EC[]'] * NT_cotribs[BB]['ATP']) + (dictionary['glucose']['EC[]'] * NT_cotribs[BB]['Glucose']) + (dictionary['glutamite']['EC[]'] * NT_cotribs[BB]['Glutamite']) + (dictionary['O2']['EC[]'] * NT_cotribs[BB]['O2']))
            elif celltype == "SC":
                lnQ -= math.log((dictionary['ATP']['SC[]'] * NT_cotribs[BB]['ATP']) + (dictionary['glucose']['SC[]'] * NT_cotribs[BB]['Glucose']) + (dictionary['glutamite']['SC[]'] * NT_cotribs[BB]['Glutamite']) + (dictionary['O2']['SC[]'] * NT_cotribs[BB]['O2']))
            elif celltype == "M":
                lnQ -= math.log((dictionary['ATP']['M[]'] * NT_cotribs[BB]['ATP']) + (dictionary['glucose']['M[]'] * NT_cotribs[BB]['Glucose']) + (dictionary['glutamite']['M[]'] * NT_cotribs[BB]['Glutamite']) + (dictionary['O2']['M[]'] * NT_cotribs[BB]['O2']))
    return lnQ


'''Omics energy'''
def get_omic_energy(all_molar_Gibbs, av_U_size):
    av_dG = statistics.mean(all_molar_Gibbs) #take the average, which is in kJ/mol
    dG_KJg = av_dG / (1 * av_U_size) # convert to kJ / dry g
    return dG_KJg


'''Mean length of the biomolecule'''
def get_mean_length (tot_bb, N_biomolecules):
    #Total number of nucleotides or amino acids divided by the number of genes/transcripts/proteins
    mean_lenght = tot_bb / N_biomolecules
    return mean_lenght



def get_Metabolites_dict(T):
    _dAS = BioMolecule('Deoxyadenosine(aq)', 251.24192, 13, T=T)
    _A = BioMolecule('Adenine(aq)', 135, 5, T=T)
    _dCS = BioMolecule('Deoxycytidine(aq)', 227.21722, 13, T=T)
    _C = BioMolecule('Cytosine(aq)', 111.102, 5, T=T)
    _dGS = BioMolecule('Deoxyguanosine(aq)', 267.24132, 13, T=T)
    _G = BioMolecule('Guanine(aq)', 151.1261, 5, T=T)
    _dTS = BioMolecule('Deoxythymidine(aq)', 242.22856, 14, T=T)
    _T = BioMolecule('Thymine(aq)', 126.11334, 6, T=T)

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
    ester_b1 = _P2O7.stdbio_formation_gibbs - (( 2 *_HPO4.stdbio_formation_gibbs) - _Water.stdbio_formation_gibbs)
    ester_b2 = _R5P.stdbio_formation_gibbs  - ((_rib.stdbio_formation_gibbs + Phosphate) - _Water.stdbio_formation_gibbs)
    ester_bond = ((ester_b1 + ester_b2) / 2)
    A_glycosidic_bond = _dAS.stdbio_formation_gibbs - ((_A.stdbio_formation_gibbs + _drib.stdbio_formation_gibbs) - _Water.stdbio_formation_gibbs)
    C_glycosidic_bond = _dCS.stdbio_formation_gibbs - ((_C.stdbio_formation_gibbs + _drib.stdbio_formation_gibbs) - _Water.stdbio_formation_gibbs)
    G_glycosidic_bond = _dGS.stdbio_formation_gibbs - ((_G.stdbio_formation_gibbs + _drib.stdbio_formation_gibbs) - _Water.stdbio_formation_gibbs)
    T_glycosidic_bond = _dTS.stdbio_formation_gibbs - ((_T.stdbio_formation_gibbs + _drib.stdbio_formation_gibbs) - _Water.stdbio_formation_gibbs)
    all_gb = [A_glycosidic_bond, C_glycosidic_bond, G_glycosidic_bond, T_glycosidic_bond]
    N_glycosidic_bond = mean(all_gb)
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



    #The dictionaries for metabolites:
    #                       0= Molecular formula
    #                       1= Molecular weight (Da)
    #                       2= dGf
    #                       3= E. coli - absolute intracellular concentration (mol/L) (Bennett et al. 2009 & Park et al. 2016)
    #                       4= Yeast - absolute intracellular concentration (mol/L) (Bennett et al. 2009 & Park et al. 2016)
    #                       5= Mammalian - absolute intracellular concentration (mol/L) (Bennett et al. 2009 & Park et al. 2016)
#
    na_Metabolites = {
    'glucose': {'name': 'C6H12O6','MW': _C6H12O6.Mr, 'dGf': _C6H12O6.stdbio_formation_gibbs, 'EC[]': 7.88E-3,'SC[]': 5.31E-03,'M[]': 6.75E-04},
    'ATP': {'name': 'C10H16N5O13P3','MW': _C10H16N5O13P3.Mr, 'dGf': _C10H16N5O13P3.stdbio_formation_gibbs, 'EC[]': 9.63E-03,'SC[]': 1.93E-03,'M[]': 4.67E-03},
    'ADP': {'name': 'C10H15N5O10P2','MW': _C10H15N5O10P2.Mr, 'dGf': _C10H15N5O10P2.stdbio_formation_gibbs, 'EC[]': 5.55E-04,'SC[]': 4.88E-04,'M[]': 5.69E-04},
    'glutamite': {'name': 'C5H10N2O3','MW': _C5H10N2O3.Mr, 'dGf': _C5H10N2O3.stdbio_formation_gibbs, 'EC[]': 3.81E-03,'SC[]': 3.55E-02,'M[]': 1.62E-02},
    'serine': {'name': 'C3H7NO3','MW': _C3H7NO3.Mr, 'dGf': _C3H7NO3.stdbio_formation_gibbs, 'EC[]': 6.8E-5,'SC[]': 6.8E-5,'M[]': 1.81E-5},
    'O2': {'name': 'O2','MW': _O2.Mr, 'dGf': _O2.stdbio_formation_gibbs, 'EC[]': '?','SC[]':'?','M[]':'?'},
    'CO2': {'name': 'CO2','MW': _CO2.Mr, 'dGf': _CO2.stdbio_formation_gibbs, 'EC[]': 7.52E-05,'SC[]': 8.16E-05,'M[]': 7.63E-03},
    'H2O': {'name': 'H2O','MW': _Water.Mr, 'dGf': _Water.stdbio_formation_gibbs, 'EC[]': 5.56E-14,'SC[]': 6.11E-12,'M[]': 1.668E-10},
    'PO4': {'name': 'PO4','MW': _PO4.Mr, 'dGf': _PO4.stdbio_formation_gibbs, 'EC[]': 2.39E-02,'SC[]': 4.93E-02,'M[]': 5.83E-03},
    'H2PO4': {'name': 'H2PO4','MW': _H2PO4.Mr, 'dGf': _H2PO4.stdbio_formation_gibbs, 'EC[]': 2.39E-02,'SC[]': 4.93E-02,'M[]': 5.83E-03},
    'phosphate': {'name': 'HPO4','MW': _HPO4.Mr, 'dGf': _HPO4.stdbio_formation_gibbs, 'EC[]': 2.39E-02,'SC[]': 4.93E-02,'M[]': 5.83E-03},
    'malonate': {'name': 'C3H2O4','MW': _C3H2O4.Mr, 'dGf': _C3H2O4.stdbio_formation_gibbs, 'EC[]': '?','SC[]': '?','M[]': '?'},
    'R5P': {'name': 'C5H11O8P','MW': _R5P.Mr, 'dGf': _R5P.stdbio_formation_gibbs, 'EC[]': 1.3E-3,'SC[]': 1.52E-04,'M[]': 7.83E-05},
    'deoxyribose': {'name': 'C5H10O4','MW': _drib.Mr, 'dGf': _drib.stdbio_formation_gibbs, 'EC[]': 3.03E-4,'SC[]': '?','M[]': '?'},
    'ribose': {'name': 'C5H10O5','MW': _rib.Mr, 'dGf': _rib.stdbio_formation_gibbs, 'EC[]': 7.83E-5,'SC[]': 1.52E-4,'M[]': '?'},
    'ester_bond': {'name': 'RCO2R','MW': 'none', 'dGf': ester_bond, 'EC[]': 'none','SC[]': 'none','M[]': 'none'},
    'glycosidic_bond': {'name': 'ROH','MW': 'none', 'dGf': N_glycosidic_bond, 'EC[]': 'none','SC[]': 'none','M[]': 'none'},
    'DNA_PBB': {'name': 'Phosphate + ester bond + Deoxyribose','MW': 'none', 'dGf': DNA_phosphate_BB, 'EC[]': 'none','SC[]': 'none','M[]': 'none'},
    'RNA_PBB': {'name': 'Phosphate + ester bond + Ribose','MW': 'none', 'dGf': RNA_phosphate_BB, 'EC[]': 'none','SC[]': 'none','M[]': 'none'},
    'phosphoester': {'name': 'Phosphate + ester bond','MW': 'none', 'dGf': phosphoester, 'EC[]': 'none','SC[]': 'none','M[]': 'none'},
    'OH': {'name': 'OH','MW': 'none', 'dGf': _OH, 'EC[]': 'none','SC[]': 'none','M[]': 'none'}}


    return na_Metabolites
