'''Last edit 27.02.2022'''

"""
Synpro, a class for calculating the energetic costs assoicated with
synthesising a proteome.

An instance describes the proteome of a cell as a set and of its
individual parts at selected environmental conditions. At present, the
only condition which varies behaviour is temperature.

"""

import sys,os, math, statistics, re
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.Seq import transcribe
from BioMolecule import BioMolecule
from copy import deepcopy
from collections import Counter

from GCAtools import *

class Synpro:
    """
    Class for calculating the energetic costs assoicated with
    synthesising a proteome.

    Properties
    ----------
    T : float
        temperature governing themodynamic constants.
    Amino_Acids : dict
        dictionary of thermodynamic data for amino acids.
    Metabolites : dict
        dictionary of thermodynamic data for metabolites.
    cell : Syncell.Cell
        instance of Cell, collecting microbial properties.
    protein_list : list
        list of protein sequences that make up the proteome.
    """

    def __init__(self, T, cell):
        self.T = T
        self.Amino_Acids = Synpro.get_Amino_Acids_dict(T)
        self.Metabolites = get_Metabolites_dict(T)
        self.cell = cell
        self.protein_list = self.read_proteome()


    def update_T(self, T):
        """
        Update the properties of this Synpro object relating to new
        temperature T.
        """
        self.T = T
        self.Amino_Acids = Synpro.get_Amino_Acids_dict(T)
        # self.Metabolites = get_Metabolites_dict(T)



    #--------------------------------
    @staticmethod
    def get_Amino_Acids_dict(T):
        """
        Return a ditionary of the twenty core amino acids, along with
        relevant 'backbone' species.

        Glycine is considered a 'Special case'.
        Each dictionary key corresponds to a list, whose indicies represent:
            0 = AA 3 letter code
            1 = Molecular weight of the AA (Da)
            2 = Chemical standard R group dGf (KJ mol)
            3 = Chemical standard AA dGf (KJ mol)
            4 = Biological standard R group dGf (KJ mol)
            5 = Biological standard AA dGf (KJ mol)
            6 = E coli absolute intraself.cellular concentration (mol/L)
                - Bennett B. et al 2009 & Park J. et al 2016
            7 = S cerevisiae absolute intraself.cellular concentration (mol/L)
                - Park J. et al 2016
            8 = Average mammalian self.cell absolute intraself.cellular concentration (mol/L)
                - Park J. et al 2016

        """

        _AABB = BioMolecule('AABB', 74.05866, 4, T=T)
        _PBB = BioMolecule('PBB', 56.04918, 2, T=T)
        _ALA = BioMolecule('Alanine(aq)', 89, 7, T=T)
        _ARG = BioMolecule('Arginine(aq)', 174, 14, T=T)
        _ASN = BioMolecule('Asparagine(aq)', 132, 8, T=T)
        _ASP = BioMolecule('Aspartic-Acid(aq)', 133, 7, T=T)
        _CYS = BioMolecule('Cysteine(aq)', 121, 7, T=T)
        _GLN = BioMolecule('Glutamic-Acid(aq)', 148, 9, T=T)
        _GLU = BioMolecule('Glutamine(aq)', 147, 10, T=T)
        _HIS = BioMolecule('Histidine(aq)', 155, 9, T=T)
        _ILE = BioMolecule('Isoleucine(aq)', 131, 13, T=T)
        _LEU = BioMolecule('Leucine(aq)', 131, 13, T=T)
        _LYS = BioMolecule('Lysine(aq)', 146, 14, T=T)
        _MET = BioMolecule('Methionine(aq)', 149, 11, T=T)
        _PHE = BioMolecule('Phenylalanine(aq)', 165, 11, T=T)
        _PRO = BioMolecule('Proline(aq)', 115, 9, T=T)
        _SER = BioMolecule('Serine(aq)', 105, 7, T=T)
        _THR = BioMolecule('Threonine(aq)', 119, 9, T=T)
        _TRP = BioMolecule('Tryptophan(aq)', 181, 12, T=T)
        _TYR = BioMolecule('Tyrosine(aq)', 181, 11, T=T)
        _VAL = BioMolecule('Valine(aq)', 117, 11, T=T)
        _GLY = BioMolecule('Glycine(aq)', 75, 5, T=T)
        _GLYlink = BioMolecule('GLY', 57, 3, T=T)
        _Water = BioMolecule('H2O(l)', 18, 2, T=T)

        AAdict = {
        'A': { 'name': 'ALA', 'MW': _ALA.Mr, 'dGf_R': _ALA.std_formation_R, 'dGf_AA': _ALA.std_formation_gibbs, 'bio_dGf_R': _ALA.stdbio_formation_R, 'bio_dGf_AA': _ALA.stdbio_formation_gibbs, 'EC[]': 2.81E-04, 'SC[]': 8.12E-05, 'M[]': 4.23E-05},
        'R': { 'name': 'ARG', 'MW': _ARG.Mr, 'dGf_R': _ARG.std_formation_R, 'dGf_AA': _ARG.std_formation_gibbs, 'bio_dGf_R': _ARG.stdbio_formation_R, 'bio_dGf_AA': _ARG.stdbio_formation_gibbs, 'EC[]': 5.69E-04, 'SC[]': 2.18E-02, 'M[]': 2.55E-04},
        'N': { 'name': 'ASN', 'MW': _ASN.Mr, 'dGf_R': _ASN.std_formation_R, 'dGf_AA': _ASN.std_formation_gibbs, 'bio_dGf_R': _ASN.stdbio_formation_R, 'bio_dGf_AA': _ASN.stdbio_formation_gibbs, 'EC[]': 5.11E-04, 'SC[]': 5.69E-03, 'M[]': 1.56E-04},
        'D': { 'name': 'ASP', 'MW': _ASP.Mr, 'dGf_R': _ASP.std_formation_R, 'dGf_AA': _ASP.std_formation_gibbs, 'bio_dGf_R': _ASP.stdbio_formation_R, 'bio_dGf_AA': _ASP.stdbio_formation_gibbs, 'EC[]': 4.23E-03, 'SC[]': 6.29E-03, 'M[]': 1.49E-02},
        'C': { 'name': 'CYS', 'MW': _CYS.Mr, 'dGf_R': _CYS.std_formation_R, 'dGf_AA': _CYS.std_formation_gibbs, 'bio_dGf_R': _CYS.stdbio_formation_R, 'bio_dGf_AA': _CYS.stdbio_formation_gibbs, 'EC[]': 8.40E-05, 'SC[]': 8.40E-05, 'M[]': 8.40E-05},
        'Q': { 'name': 'GLN', 'MW': _GLN.Mr, 'dGf_R': _GLN.std_formation_R, 'dGf_AA': _GLN.std_formation_gibbs, 'bio_dGf_R': _GLN.stdbio_formation_R, 'bio_dGf_AA': _GLN.stdbio_formation_gibbs, 'EC[]': 3.81E-03, 'SC[]': 3.55E-02, 'M[]': 1.62E-02},
        'E': { 'name': 'GLU', 'MW': _GLU.Mr, 'dGf_R': _GLU.std_formation_R, 'dGf_AA': _GLU.std_formation_gibbs, 'bio_dGf_R': _GLU.stdbio_formation_R, 'bio_dGf_AA': _GLU.stdbio_formation_gibbs, 'EC[]': 9.60E-02, 'SC[]': 3.91E-02, 'M[]': 4.36E-02},
        'H': { 'name': 'HIS', 'MW': _HIS.Mr, 'dGf_R': _HIS.std_formation_R, 'dGf_AA': _HIS.std_formation_gibbs, 'bio_dGf_R': _HIS.stdbio_formation_R, 'bio_dGf_AA': _HIS.stdbio_formation_gibbs, 'EC[]': 6.76E-05, 'SC[]': 6.76E-05, 'M[]': 4.10E-04},
        'I': { 'name': 'ILE', 'MW': _ILE.Mr, 'dGf_R': _ILE.std_formation_R, 'dGf_AA': _ILE.std_formation_gibbs, 'bio_dGf_R': _ILE.stdbio_formation_R, 'bio_dGf_AA': _ILE.stdbio_formation_gibbs, 'EC[]': 1.52E-04, 'SC[]': 3.53E-04, 'M[]': 1.66E-03},
        'L': { 'name': 'LEU', 'MW': _LEU.Mr, 'dGf_R': _LEU.std_formation_R, 'dGf_AA': _LEU.std_formation_gibbs, 'bio_dGf_R': _LEU.stdbio_formation_R, 'bio_dGf_AA': _LEU.stdbio_formation_gibbs, 'EC[]': 1.52E-04, 'SC[]': 3.53E-04, 'M[]': 1.66E-03},
        'K': { 'name': 'LYS', 'MW': _LYS.Mr, 'dGf_R': _LYS.std_formation_R, 'dGf_AA': _LYS.std_formation_gibbs, 'bio_dGf_R': _LYS.stdbio_formation_R, 'bio_dGf_AA': _LYS.stdbio_formation_gibbs, 'EC[]': 4.05E-04, 'SC[]': 5.16E-03, 'M[]': 5.06E-04},
        'M': { 'name': 'MET', 'MW': _MET.Mr, 'dGf_R': _MET.std_formation_R, 'dGf_AA': _MET.std_formation_gibbs, 'bio_dGf_R': _MET.stdbio_formation_R, 'bio_dGf_AA': _MET.stdbio_formation_gibbs, 'EC[]': 1.45E-04, 'SC[]': 1.91E-04, 'M[]': 6.19E-04},
        'F': { 'name': 'PHE', 'MW': _PHE.Mr, 'dGf_R': _PHE.std_formation_R, 'dGf_AA': _PHE.std_formation_gibbs, 'bio_dGf_R': _PHE.stdbio_formation_R, 'bio_dGf_AA': _PHE.stdbio_formation_gibbs, 'EC[]': 1.82E-05, 'SC[]': 2.73E-04, 'M[]': 7.97E-04},
        'P': { 'name': 'PRO', 'MW': _PRO.Mr, 'dGf_R': _PRO.std_formation_R, 'dGf_AA': _PRO.std_formation_gibbs, 'bio_dGf_R': _PRO.stdbio_formation_R, 'bio_dGf_AA': _PRO.stdbio_formation_gibbs, 'EC[]': 3.85E-04, 'SC[]': 1.36E-03, 'M[]': 1.23E-03},
        'S': { 'name': 'SER', 'MW': _SER.Mr, 'dGf_R': _SER.std_formation_R, 'dGf_AA': _SER.std_formation_gibbs, 'bio_dGf_R': _SER.stdbio_formation_R, 'bio_dGf_AA': _SER.stdbio_formation_gibbs, 'EC[]': 1.13E-03, 'SC[]': 3.87E-03, 'M[]': 4.86E-03},
        'T': { 'name': 'THR', 'MW': _THR.Mr, 'dGf_R': _THR.std_formation_R, 'dGf_AA': _THR.std_formation_gibbs, 'bio_dGf_R': _THR.stdbio_formation_R, 'bio_dGf_AA': _THR.stdbio_formation_gibbs, 'EC[]': 1.26E-03, 'SC[]': 6.69E-03, 'M[]': 6.69E-03},
        'W': { 'name': 'TRP', 'MW': _TRP.Mr, 'dGf_R': _TRP.std_formation_R, 'dGf_AA': _TRP.std_formation_gibbs, 'bio_dGf_R': _TRP.stdbio_formation_R, 'bio_dGf_AA': _TRP.stdbio_formation_gibbs, 'EC[]': 1.21E-05, 'SC[]': 5.55E-05, 'M[]': 1.80E-04},
        'Y': { 'name': 'TYR', 'MW': _TYR.Mr, 'dGf_R': _TYR.std_formation_R, 'dGf_AA': _TYR.std_formation_gibbs, 'bio_dGf_R': _TYR.stdbio_formation_R, 'bio_dGf_AA': _TYR.stdbio_formation_gibbs, 'EC[]': 2.89E-05, 'SC[]': 2.48E-04, 'M[]': 8.88E-04},
        'V': { 'name': 'VAL', 'MW': _VAL.Mr, 'dGf_R': _VAL.std_formation_R, 'dGf_AA': _VAL.std_formation_gibbs, 'bio_dGf_R': _VAL.stdbio_formation_R, 'bio_dGf_AA': _VAL.stdbio_formation_gibbs, 'EC[]': 4.02E-03, 'SC[]': 2.50E-03, 'M[]': 1.44E-03},
        'G': { 'name': 'GLY', 'MW': _GLY.Mr, 'dGf_R': _GLY.std_formation_R, 'dGf_AA': _GLY.std_formation_gibbs, 'bio_dGf_R': _GLY.stdbio_formation_R, 'bio_dGf_AA': _GLY.stdbio_formation_gibbs, 'EC[]': 3.71E-03, 'SC[]': 3.71E-03, 'M[]': 3.71E-03},
        '_G': { 'name': '-GLY-', 'MW': _GLYlink.Mr, 'dGf_AA': _GLYlink.std_formation_gibbs, 'bio_dGf_AA': _GLYlink.stdbio_formation_gibbs, 'EC[]': '?', 'SC[]': '?', 'M[]':'?'},
        'AABB': { 'name': 'AABB', 'MW': _AABB.Mr, 'dGf_AA': _AABB.std_formation_gibbs, 'bio_dGf_AA': _AABB.stdbio_formation_gibbs, 'EC[]': '?', 'SC[]': '?', 'M[]': '?'},
        'PBB': { 'name': 'PBB', 'MW': _PBB.Mr, 'dGf_AA': _PBB.std_formation_gibbs, 'bio_dGf_AA': _PBB.stdbio_formation_gibbs, 'EC[]': '?', 'SC[]': '?', 'M[]': '?'},
        'H2O': {'name': 'H2O','MW': _Water.Mr, 'dGf_AA': _Water.std_formation_gibbs, 'EC[]': 1,'SC[]': 1,'M[]': 1, 'bio_dGf_AA': _Water.stdbio_formation_gibbs}}

        return AAdict




    '''Formation energy of the proteins'''
    @staticmethod
    def get_dGf_P (protein, AA, standard, R_standard):
        sigma_AA_R = 0
        local_AA_R_count = []
        sigma_Gly = 0
        local_Gly_count = []

        PBB = AA['PBB'][standard]
        AABB = AA['AABB'][standard]

        for aa in protein:
            if aa == 'G':
                sigma_Gly += AA['_G'][standard]
                local_Gly_count.append(sigma_Gly)
            else:
                if aa in AA.keys():
                    sigma_AA_R += AA[aa][R_standard]
                    local_AA_R_count.append(sigma_AA_R)

        Nn = len(local_AA_R_count) + len(local_Gly_count)
        NGly = len(local_Gly_count)

        #Equation 4 from self.cell Biosynthesis Energetics
        protein_dGf = (AABB + ((Nn - NGly - 1) * PBB)) + sigma_AA_R + sigma_Gly
        return protein_dGf


    '''reaction energy of the proteins'''
    def get_dGr_P(self, protein, dGf_protein, AA, standard):
        Nn = len(protein)
        H2O = self.Amino_Acids['H2O'][standard]
        sigma_AA = 0

        for aa in protein:
            if aa in AA.keys():
                sigma_AA += AA[aa][standard]

        protein_dGr = (((Nn - 1) * H2O) + dGf_protein) - sigma_AA
        return protein_dGr


    def read_proteome(self):
        """read the proteome fasta file of a self.cell object and return a list of sequences """
        protein_list = []
        '''-Read Sequence-'''
        fasta_p = open(self.cell.get_Proteome(),'r')
        for record in SeqIO.parse(fasta_p,'fasta'):
            seqs = str(record.seq)
            protein_list.append(str(re.sub('[^a-zA-Z]','',seqs)))

        return protein_list



    def proteome_synthesis(self):
        """
        Calculate and return the free energy of synthesising one dry g of the
        proteome of this instance's Cell.

        Two energies are returned, the first using chemical standard free
        energies and the second using biological standard free energies.
        """

        #--- Number of amino acids on each protein and total weight---
        '''-----Proteins-----'''
        tot_AA = 0
        proteome_weight = 0
        #---Energy of each strand (two conditions )--
        # Chemical standard
        protein_dGf_list = []
        protein_dGr_list = []
        protein_reactionGs = []
        # Biological standard
        protein_bio_dGf_list = []
        protein_bio_dGr_list = []
        protein_bio_reactionGs = []
        #---- Amino acids count
        for protein in self.protein_list:
            '''General characteristics of the gene'''
            # Amino acid count and protein lenght
            cur_dict = deepcopy(self.Amino_Acids)
            AA_count = Counter(protein)
            for aa, dictionary in cur_dict.items():
                dictionary['count'] = AA_count.get(aa, 0)

            NAA = len(protein) # Number of amino acids in each protein
            N_links = len(protein) - 1 # Number of peptide links in each protein
            NGly = cur_dict['G']['count'] # Number of glycines
            tot_AA += NAA # Total number of amino acids in the proteome
            protein_meanlength = get_mean_length(tot_AA, len(self.protein_list)) # Average protein length in this proteome

            # Proteins and proteome weight
            protein_weight = get_weight(cur_dict) - (N_links * self.Metabolites['H2O']['MW']) # Weight of each protein
            proteome_weight += protein_weight # Weight of the proteome considering just one copy of each sequence


            '''---------------------'''
            '''Protein Formation energy'''
            # Chemical standard
            dGf_protein = Synpro.get_dGf_P(protein, cur_dict, 'dGf_AA', 'dGf_R')
            protein_dGf_list.append(dGf_protein)
            #print('Formation energy of the protein at chemical standard', dGf_protein, 'at', T, 'K')
            # Biological standard
            bio_dGf_protein = Synpro.get_dGf_P(protein, cur_dict, 'bio_dGf_AA', 'bio_dGf_R')
            protein_bio_dGf_list.append(bio_dGf_protein)
            #print('Formation energy of the protein at biological standard', dGf_protein, 'at', T, 'K')

            '''---------------------'''
            '''Protein Reaction energy'''
            # Chemical standard
            dGr_protein = self.get_dGr_P(protein, dGf_protein, cur_dict, 'dGf_AA')
            protein_dGr_list.append(dGr_protein)
            # print('Reaction energy of the protein at chemical standard', dGr_protein, 'at', T, 'K')
            # Biological standard
            bio_dGr_protein = self.get_dGr_P(protein, bio_dGf_protein, cur_dict, 'bio_dGf_AA',)
            protein_bio_dGr_list.append(bio_dGr_protein)
            # print('Reaction energy of the protein at biological standard', dGr_protein, 'at', T, 'K')

        '''Molarity'''
        av_aa_mw = proteome_weight / tot_AA # Average amino acid  weight (Da)
        av_protein_w = proteome_weight / len(self.protein_list) # Average protein weight (Da)
        mean_prot_length = av_protein_w / av_aa_mw
        av_aa_mass = av_aa_mw / 6.02E+23 #Average amino acid weight divided by avogadro's number
        aa_moles = self.cell.get_proteins_weight() / av_aa_mass
        protein_moles = aa_moles / mean_prot_length
        protein_molarity = protein_moles * self.cell.get_volume()



        '''---------------------'''
        '''Reaction quotient'''
        R = 0.008314472 #Gas constant KJ mol K
        for protein, protein_dGr, protein_bio_dGr in zip(self.protein_list, protein_dGr_list, protein_bio_dGr_list):
            protein_lnQ = get_lnQ(protein_molarity, protein, len(self.protein_list), cur_dict, self.cell.get_model())

            '''Molar Gibbs Energy of each protein'''
            '''Chemical standard'''
            protein_dG = (protein_dGr + (R * self.T * protein_lnQ))
            protein_reactionGs.append(protein_dG)

            '''Biological standard'''
            protein_bio_dG = (protein_bio_dGr + (R * self.T * protein_lnQ))
            protein_bio_reactionGs.append(protein_bio_dG)


        '''Proteome energy (KJ/g)'''
        '''Chemical standard'''
        proteome_energy = get_omic_energy(protein_reactionGs, av_protein_w)

        # print('proteome energy chemical standard', proteome_energy, 'KJ/g at ', T, 'K')
        '''Biological standard'''
        proteome_bio_energy = get_omic_energy(protein_bio_reactionGs, av_protein_w)

        # print('proteome energy biological standard', proteome_bio_energy, 'KJ/g at ', T, 'K')

        return proteome_energy, proteome_bio_energy
