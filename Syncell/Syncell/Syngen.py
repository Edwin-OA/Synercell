'''Last edit 03.03.2022 - Edwin'''

"""
Syngen, a class for calculating the energetic costs assoicated with
synthesising a genome and transcriptome.

An instance describes the genome of a cell as a set and of its
individual parts at selected environmental conditions. At present, the
only condition which varies behaviour is temperature.

"""

import sys,os, math, statistics, re
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.Seq import transcribe
from BioMolecule import BioMolecule

import time
from copy import deepcopy
from collections import Counter

from GCAtools import *

class Syngen:
    """
    Class for calculating the energetic costs assoicated with
    synthesising the genome and transcriptome of a cell.

    Properties
    ----------
    T : float
        temperature governing themodynamic constants.
    Nucleic_Acids : dict
        dictionary of thermodynamic data for nucleic acids.
    Metabolites : dict
        dictionary of thermodynamic data for metabolites.
    cell : Syncell.Cell
        instance of Cell, collecting microbial properties.
    gene_list : list
        list of DNA sequences that make up the proteome.
    mRNA_list : list
        list of mRNA sequences that make up the transcriptome.
    """

    def __init__(self, T, cell):
        self.T = T
        self.Nucleic_Acids = Syngen.get_Nucleic_Acids_dict(T)
        self.Metabolites = get_Metabolites_dict(T)
        self.cell = cell

        gene_list, mRNA_list = self.read_genome()
        self.gene_list = gene_list
        self.mRNA_list = mRNA_list

    def update_T(self, T):
        """
        Update the properties of this Syngen object relating to new
        temperature T.
        """
        self.T = T
        self.Nucleic_Acids = Syngen.get_Nucleic_Acids_dict(T)
        self.Metabolites = get_Metabolites_dict(T)

    #--------------------------------
    #Syngen dictionaries
    @staticmethod
    def get_Nucleic_Acids_dict(T):
        """
        This method returns a nested dictionary in the format:
            {'DNA':dict1, 'RNA':dict2}

        where dict1 and dict2 consist of 'A', 'C', 'G', 'T' and 'U'
        (as appropriate) BioMolecules and their properties at temperature
        T in Kelvin.

        The dictionaries contain lists with the following incidies:
                0 = Nucleotide name
                1  Nucleotide's molecular weight (Da)
        Formation energies -
                2 = dGf - base (database)
                3 = dGf - Nucleoside (database)
                4 = dGf - Nucleotide (database)
                5 = dGf - Nucleotide ion (database)
        Nucelosides intracellular concentrations -
                6 = E. coli  (mol/L) (Bennett et al. 2009 & Park et al. 2016)
                7 = Yeast  (mol/L) (Bennett et al. 2009 & Park et al. 2016)
                8 = Mammalian (mol/L) (Bennett et al. 2009 & Park et al. 2016)
        Nucleotides intracellular concentrations -
                9 = E. coli  (mol/L) (Bennett et al. 2009 & Park et al. 2016)
                10 = Yeast  (mol/L) (Bennett et al. 2009 & Park et al. 2016)
                11 = Mammalian (mol/L) (Bennett et al. 2009 & Park et al. 2016)

        """

        '''A'''
        # BioMolecules are declared with (name, Mr, no. H, T, charge)
        #Base
        _A = BioMolecule('Adenine(aq)', 135, 5, T=T)
        #DNA Nucleoside
        _dAS = BioMolecule('Deoxyadenosine(aq)', 251.24192, 13, T=T)
        #RNA nucleoside
        _AS = BioMolecule('Adenosine(aq)', 267.24132, 13, T=T)
        #DNA Nucleotide
        _dAMP2 = BioMolecule('d+H2AMP-(aq)', 300.24806, 14, T=T)
        _dAMP = BioMolecule('dHAMP-', 330.213881, 13, T=T, z=-1)
        #ion
        _dAMPion = BioMolecule('dAMP2-', 329.205941, 12, T=T, z=-2)
        #RNA nucleotide
        #_AMP = BioMolecule('+H2AMP-(aq)', 347.221221, 14, T=T)
        _AMP = BioMolecule('HAMP-', 346.205341, 13, T=T, z=-1)
        #ion
        _AMPion = BioMolecule('AMP2-', 345.205341, 12, T=T, z=-2)
        '''C'''
        # BioMolecules are declared with (name, Mr, no. H, T, charge)
        #Base
        _C = BioMolecule('Cytosine(aq)', 111.102, 5, T=T)
        #DNA Nucleoside
        _dCS = BioMolecule('Deoxycytidine(aq)', 227.21722, 13, T=T)
        #RNA nucleoside
        _CS = BioMolecule('Cytidine(aq)', 243.21662, 13, T=T)
        #DNA Nucleotide
        #_dCMP = BioMolecule('d+H2CMP-(aq)', 307.197121, 14, T=T)
        _dCMP = BioMolecule('dHCMP-', 306.189181, 13, T=T, z=-1)
        #ion
        _dCMPion = BioMolecule('dCMP2-', 305.181241, 12, T=T, z=-2)
        #RNA nucleotide
        #_CMP = BioMolecule('+H2CMP-(aq)', 323.196521, 14, T=T)
        _CMP = BioMolecule('HCMP-', 322.188581, 13, T=T, z=-1)
        #ion
        _CMPion = BioMolecule('CMP2-', 321.180641, 12, T=T, z=-2)

        '''G'''
        # BioMolecules are declared with (name, Mr, no. H, T, charge)
        #Base
        _G = BioMolecule('Guanine(aq)', 151.1261, 5, T=T)
        #DNA nucleoside
        _dGS = BioMolecule('Deoxyguanosine(aq)', 267.24132, 13, T=T)
        #RNA nucleoside
        _GS = BioMolecule('Guanosine(aq)', 283.24072, 13, T=T)
        #DNA nucleotide
        #_dGMP = BioMolecule('d+H2GMP-(aq)', 316.24746, 14, T=T)
        _dGMP = BioMolecule('dHGMP-', 346.213281, 13, T=T, z=-1)
        #ion
        _dGMPion = BioMolecule('dGMP2-', 345.205341, 12, T=T, z=-2)
        #RNA nucleotide
        #_GMP = BioMolecule('+H2GMP-(aq)', 363.220621, 14, T=T)
        _GMP = BioMolecule('HGMP-', 362.212681, 13, T=T, z=-1)
        #ion
        _GMPion = BioMolecule('GMP2-', 361.204741, 12, T=T, z=-2)

        '''T'''
        # BioMolecules are declared with (name, Mr, no. H, T, charge)
        #Base
        _T = BioMolecule('Thymine(aq)', 126.11334, 6, T=T)
        #DNA nucleoside
        _dTS = BioMolecule('Deoxythymidine(aq)', 242.22856, 14, T=T)
        #DNA nucleotide
        #_dTMP = BioMolecule('d+H2TMP-(aq)', 322.208461, 15, T=T)
        _dTMP =BioMolecule('dHTMP-', 321.200521, 14, T=T, z=-1)
        #ion
        _dTMPion =BioMolecule('dTMP2-', 320.192581, 13, T=T, z=-2)

        '''U'''
        # BioMolecules are declared with (name, Mr, no. H, T, charge)
        #Base
        _U = BioMolecule('Uracil(aq)', 112.086, 4, T=T)
        #RNA nucleoside
        _US = BioMolecule('Uridine(aq)', 244.20138, 12, T=T)
        #RNA nucleotide
        #_UMP = BioMolecule('+H2UMP-(aq)', 324.181281, 13, T=T)
        _UMP = BioMolecule('HUMP-', 323.173341, 12, T=T, z=-1)
        #ion
        _UMPion = BioMolecule('UMP2-', 322.165401, 11, T=T, z=-2)

        '''Metabolites'''
        _H2O = BioMolecule('H2O(l)', 18, 2, T=T)
        H2O = _H2O.stdbio_formation_gibbs
        _ADP = BioMolecule('+H3ADP1-(aq)', 427.201122, 15, T=T)
        ADP = _ADP.stdbio_formation_gibbs
        _ATP = BioMolecule('+H4ATP-(aq)', 507.181023, 16, T=T)
        ATP = _ATP.stdbio_formation_gibbs
        _CO2 = BioMolecule('CO2(aq)', 44.0095, 0, T=T)
        CO2 = _CO2.stdbio_formation_gibbs
        _O2 = BioMolecule('O2(aq)', 31.9988, 0, T=T)
        O2 = _O2.stdbio_formation_gibbs
        _Glucose = BioMolecule('Glucose(aq)', 180.15588, 12, T=T)
        Glucose = _Glucose.stdbio_formation_gibbs
        _Glutamine = BioMolecule('Glutamine(aq)', 146.1445, 10, T=T)
        Glutamine = _Glutamine.stdbio_formation_gibbs


        DNAdict = {
          'A': {'name': 'dAMP', 'MW': _dAMP.Mr,
          'dGf[base]': _A.stdbio_formation_gibbs, 'dGf[nucleoside]': _dAS.stdbio_formation_gibbs, 'dGf[nucleotide]': _dAMP.stdbio_formation_gibbs, 'dGf[nucleotide_i]': _dAMPion.stdbio_formation_gibbs,
          'EC[NS]': 2.82E-06, 'SC[NS]': 2.82E-06, 'M[NS]': 2.82E-06,
          'EC[]': 8.84e-6, 'SC[]': 8.84e-6, 'M[]': 1.68e-5,
          'Glucose[dGf]': (Glucose * 0.5), 'ATP[dGf]': (ATP * 1), 'Glutamine[dGf]': (Glutamine * 2.5), 'O2[dGf]': (O2 * 6.5), 'ADP[dGf]': (ADP * 1), 'H2O[dGf]': (H2O * 9.5), 'CO2[dGf]': (CO2 * 5.5)},

          'C': {'name': 'dCMP', 'MW': _dCMP.Mr,
          'dGf[base]': _C.stdbio_formation_gibbs, 'dGf[nucleoside]': _dCS.stdbio_formation_gibbs, 'dGf[nucleotide]': _dCMP.stdbio_formation_gibbs, 'dGf[nucleotide_i]': _dCMPion.stdbio_formation_gibbs,
          'EC[NS]': 2.82E-06, 'SC[NS]': 2.82E-06, 'M[NS]': 2.82E-06,
          'EC[]': 3.71e-05, 'SC[]': 3.71e-05, 'M[]': 3.71e-05,
          'Glucose[dGf]': (Glucose * 0.5), 'ATP[dGf]': (ATP * 1), 'Glutamine[dGf]': (Glutamine * 1.5), 'O2[dGf]': (O2 * 2), 'ADP[dGf]': (ADP * 1), 'H2O[dGf]': (H2O * 4.5), 'CO2[dGf]': (CO2 * 1.5)},

          'G': {'name': 'dGMP', 'MW': _dGMP.Mr,
          'dGf[base]': _G.stdbio_formation_gibbs, 'dGf[nucleoside]': _dGS.stdbio_formation_gibbs, 'dGf[nucleotide]': _dGMP.stdbio_formation_gibbs, 'dGf[nucleotide_i]': _dGMPion.stdbio_formation_gibbs,
          'EC[NS]': 5.22E-07,'SC[NS]': 5.22E-07, 'M[NS]': 5.22E-07,
          'EC[]': 5.1e-5, 'SC[]': 5.1e-5,'M[]': 5.1e-5,
          'Glucose[dGf]': (Glucose * 0.5), 'ATP[dGf]': (ATP * 1), 'Glutamine[dGf]': (Glutamine * 2.5), 'O2[dGf]': (O2 * 7), 'ADP[dGf]': (ADP * 1), 'H2O[dGf]': (H2O * 9.5), 'CO2[dGf]': (CO2 * 5.5)},

          'T': {'name': 'dTMP', 'MW': _dTMP.Mr,
          'dGf[base]': _T.stdbio_formation_gibbs, 'dGf[nucleoside]': _dTS.stdbio_formation_gibbs, 'dGf[nucleotide]': _dTMP.stdbio_formation_gibbs, 'dGf[nucleotide_i]': _dTMPion.stdbio_formation_gibbs,
          'EC[NS]': 5.22E-07,'SC[NS]': 5.22E-07, 'M[NS]': 5.22E-07,
          'EC[]': 1.18e-5, 'SC[]': 1.18e-5, 'M[]': 1.18e-5,
          'Glucose[dGf]': (Glucose * 0.5), 'ATP[dGf]': (ATP * 3), 'Glutamine[dGf]': (Glutamine * 3.5), 'O2[dGf]': (O2 * 2), 'ADP[dGf]': (ADP * 4), 'H2O[dGf]': (H2O * 7.5), 'CO2[dGf]': (CO2 * 0.5)}}


        #Dictionary for RNA
        RNAdict = {
          'A': {'name': 'AMP', 'MW': _AMP.Mr,
          'dGf[base]': _A.stdbio_formation_gibbs, 'dGf[nucleoside]': _AS.stdbio_formation_gibbs, 'dGf[nucleotide]': _AMP.stdbio_formation_gibbs, 'dGf[nucleotide_i]': _AMPion.stdbio_formation_gibbs,
          'EC[NS]': 1.31E-07, 'SC[NS]': 1.31E-07, 'M[NS]': 1.31E-07,
          'EC[]': 2.8e-4, 'SC[]': 8.12e-5, 'M[]': 4.23e-5},
          'C': {'name': 'CMP', 'MW': _CMP.Mr,
          'dGf[base]': _C.stdbio_formation_gibbs, 'dGf[nucleoside]': _CS.stdbio_formation_gibbs, 'dGf[nucleotide]': _CMP.stdbio_formation_gibbs, 'dGf[nucleotide_i]': _CMPion.stdbio_formation_gibbs,
          'EC[NS]': 1.41E-05, 'SC[NS]': 1.41E-05, 'M[NS]': 1.41E-05,
          'EC[]': 3.6e-4, 'SC[]': 5.18e-6, 'M[]': 1.18e-5},
          'G': {'name': 'GMP', 'MW': _GMP.Mr,
          'dGf[base]': _G.stdbio_formation_gibbs, 'dGf[nucleoside]': _GS.stdbio_formation_gibbs, 'dGf[nucleotide]': _GMP.stdbio_formation_gibbs, 'dGf[nucleotide_i]': _GMPion.stdbio_formation_gibbs,
          'EC[NS]': 2.09E-03, 'SC[NS]': 2.09E-03, 'M[NS]': 2.09E-03,
          'EC[]': 2.37e-5, 'SC[]': 1.02e-5, 'M[]': 1.81e-5},
          'U': {'name': 'UMP', 'MW': _UMP.Mr,
          'dGf[base]': _U.stdbio_formation_gibbs, 'dGf[nucleoside]': _US.stdbio_formation_gibbs, 'dGf[nucleotide]': _UMP.stdbio_formation_gibbs, 'dGf[nucleotide_i]': _UMPion.stdbio_formation_gibbs,
          'EC[NS]': 1.62E-06, 'SC[NS]': 1.62E-06, 'M[NS]': 1.35E-06,
          'EC[]': 1.45E-05, 'SC[]': 1.45E-05, 'M[]': 1.45E-05}}

        Nucleic_Acids = {'RNA': RNAdict, 'DNA': DNAdict}
        return Nucleic_Acids




    #--------Transcription------
    @staticmethod
    def transcription(ORF):
        mRNAlst = [] # save the fragments in a list
        mRNA= ''
        foundStart = False
        foundEnd = False
        for i in range(0, len(ORF), 3):
            codon= "".join(ORF[i:i+3]) # Define a codon
            if codon == 'ATG' and not foundStart:
                foundStart = True
                foundEnd = False # start recording until it finds TAG.
            if foundStart and not foundEnd:
                cc=transcribe(codon)
                mRNA = mRNA + cc
            if codon == 'TAA': # First possible stop codon
                foundEnd = True
            elif codon == 'TGA': # Second possible stop codon
                foundEnd = True
            elif codon == 'TAG': # Third possible stop codon
                foundEnd = True
                foundStart = False # Start looking again
                mRNAlst.append(mRNA) # save what we have to the list
                mRNA='' # reset for starting again

        # if we get to the end with no TAG, add everything saved to mRNA so far
        if mRNA != '':
            mRNAlst.append(mRNA)
        return(mRNAlst)




    '''Formation energy of the nucleic acids'''

    #For all methods
    def get_dGf_na(self, molecule, dictionary, type, method):
        #Sum of the formation energy of all the nucleotides-ion in the gene
        Nn = len(molecule)
        P_bonds = Nn - 1
        hydrogen_phosphate = self.Metabolites['phosphate']['dGf']
        ester_bond = self.Metabolites['ester_bond']['dGf']
        glyco_bond = self.Metabolites['glycosidic_bond']['dGf']
        OH = self.Metabolites['OH']['dGf']
        H2O = self.Metabolites['H2O']['dGf']
        if type == 'DNA':
            sugar = self.Metabolites['deoxyribose']['dGf']
        if type == 'RNA':
            sugar = self.Metabolites['ribose']['dGf']
        sigma = 0
        #-------------------
        if method == 'chain':
            for n in molecule:
                if n in dictionary.keys():
                    sigma += dictionary[n]['dGf[nucleotide_i]']

        if method == 'block':
            for n in molecule:
                if n in dictionary.keys():
                    sigma += (((hydrogen_phosphate + ester_bond + sugar + glyco_bond) - (2 * H2O)) + dictionary[n]['dGf[base]'])

        if method == 'block2':
            for n in molecule:
                if n in dictionary.keys():
                    sigma += (((hydrogen_phosphate + ester_bond) - (H2O)) + dictionary[n]['dGf[nucleoside]'])


        #Equation 6 from Cell Biosynthesis Energetics
        na_dGf = ((P_bonds) * (ester_bond - OH)) +  sigma
        return na_dGf


    '''Reaction energy of the nucleic acids'''

    def get_dGr_na (self, na, NA, dGf_na, method):
        #na: Nucleic Acid (molecule), NA: dictionary the na, dGf_na: formation energy of the na, method: minimum or metabolites
        #Reactants and products
        H2O = self.Metabolites['H2O']['dGf']
        Nn = len(na)
        P_bonds =  Nn - 1
        sigma_reactants = 0
        sigma_products = 0

        if method == 'minimum':
            for n in na:
                if n in NA.keys():
                    sigma_reactants += NA[n]['dGf[nucleotide]']

        #Equation 7 from Cell Biosynthesis Energetics
        na_dGr = ((P_bonds * H2O) + dGf_na) - sigma_reactants
        return na_dGr

        if method == 'metabolites':
            for n in na:
                if n in NA.keys():
                    sigma_reactants += NA[n]['Glucose[dGf]'] + NA[n]['ATP[dGf]'] + NA[n]['Glutamine[dGf]'] + NA[n]['O2[dGf]']
                    sigma_products += NA[n]['ADP[dGf]'] + NA[n]['H2O[dGf]'] + NA[n]['CO2[dGf]']

        #Equation 10 from Cell Biosynthesis Energetics
        na_dGr = ((P_bonds * H2O) + dGf_na + sigma_products) - (sigma_reactants)
        return na_dGr


    def read_genome(self):
        """read the proteome fasta file of a self.cell
        object and return two lists of sequences, the first of the genome
        and the second of the mRNA.  """

        Seqname = []
        gene_list = []
        records = []
        raw_RNAs = []
        mRNA_list = []
        noncoding_RNAs = []

        '''-Read Sequence-'''
        fasta = open(self.cell.get_Genome(),'r')
        for record in SeqIO.parse(fasta,'fasta'):
            time.sleep(0.01)
            Seqname.append(record.id) #Name of the gene
            gene_list.append(str(record.seq)) #Genes in a list of strings
            records.append((record.seq)) # Genes in a list of int to be used by biopyhton

        '''-Complement Sequence-'''
        for record in records:
            time.sleep(0.01)
            complementary=record.complement()
            gene_list.append(str(complementary))

        '''-Transcribe Sequence-'''
        for gene in records:
            time.sleep(0.01)
            RNA = Syngen.transcription(gene)
            if RNA != '':
                raw_RNAs.extend(RNA) #  Add a list on the end to be cleared
        transcripts = [x for x in raw_RNAs if x]

        '''-Clean transcripts-'''
        for rna in transcripts:
            if len(rna) > 18 and len(rna) < 750: #Shea J. Andrews and Joseph A. Rothnagel. Nature Reviews  2014. "The smallest translated coding sORF described so far is six codons long."
                mRNA_list.append(rna)
            else:
                noncoding_RNAs.append(rna)

        return gene_list, mRNA_list


    def genome_synthesis(self):
        """
        Calculate and return the free energy of synthesising one dry g of the
        genome of this instance's Cell.
        """

        #Nucleotide count for each gene
        #--- Number of nucleotides on each strand and total weight---
        '''-----DNA-----'''
        tot_n = 0
        genome_weight = 0
        #---Energy of each strand (Chain method )--
        DNA_dGf_list = []
        DNA_dGr_list = []
        DNA_reactionGs = []

        for gene in self.gene_list:

            '''General characteristics of the gene'''
            #Nucleotide count and gene lenght
            DNA_cur_dict = deepcopy(self.Nucleic_Acids['DNA'])
            DNA_n_count = Counter(gene)

            for n, dictionary in DNA_cur_dict.items():
                dictionary['count'] = DNA_n_count.get(n ,0)

            Nn = len(gene)
            P_bonds = len(gene) - 1
            tot_n += Nn
            gene_meanlength = get_mean_length(tot_n, len(self.gene_list))

            # Gene and genome weight
            gene_weight = get_weight(DNA_cur_dict) - (P_bonds * self.Metabolites['H2O']['MW'])
            genome_weight += gene_weight # Weight of the genome (Da)


            '''---------------------'''
            '''Gene Formation energy'''
            dGf_gene = self.get_dGf_na(gene, DNA_cur_dict,'DNA','chain')

            '''---------------------'''
            '''Gene Reaction energy'''
            dGr_gene = self.get_dGr_na(gene, DNA_cur_dict, dGf_gene, 'minimum')
            # dGr_gene = self.get_dGr_na(gene, DNA_cur_dict, dGf_gene, 'metabolites')

            DNA_dGr_list.append(dGr_gene)

        '''Molarity'''
        nucleotide_av_mw = genome_weight / tot_n # Average nucleotide weight (Da)
        av_gene_w = genome_weight / len(self.gene_list) # Average gene weight (Da)
        DNA_moles = self.cell.get_DNA_weight() / nucleotide_av_mw # Moles of DNA
        DNA_molarity = DNA_moles / self.cell.get_volume()

        '''Reqction quotient'''
        R = 0.008314472 #Gas constant KJ mol K
        for gene, gene_dGr in zip(self.gene_list, DNA_dGr_list):
            gene_lnQ = get_lnQ(DNA_molarity, gene, len(self.gene_list), self.Nucleic_Acids['DNA'], self.cell.get_model())
            # gene_lnQ = get_alt_lnQ(DNA_molarity, gene, len(self.gene_list), self.Metabolites, self.cell.get_model())

            '''Molar Gibbs Energy of each gene'''
            DNA_dG = (gene_dGr + (R * self.T * gene_lnQ))
            DNA_reactionGs.append(DNA_dG)

        '''Genome energy (KJ/g)'''
        genome_energy = get_omic_energy(DNA_reactionGs, av_gene_w)

        return genome_energy

    def genome_synthesis1(self):
        """
        Calculate and return the free energy of synthesising one dry g of the
        genome of this instance's Cell.
        """

        #Nucleotide count for each gene
        #--- Number of nucleotides on each strand and total weight---
        '''-----DNA-----'''
        tot_n = 0
        genome_weight = 0
        #---Energy of each strand (Chain method )--
        DNA_dGf_list = []
        DNA_dGr_list1 = []
        DNA_reactionGs1 = []

        for gene in self.gene_list:

            '''General characteristics of the gene'''
            #Nucleotide count and gene lenght
            DNA_cur_dict = deepcopy(self.Nucleic_Acids['DNA'])
            DNA_n_count = Counter(gene)

            for n, dictionary in DNA_cur_dict.items():
                dictionary['count'] = DNA_n_count.get(n ,0)

            Nn = len(gene)
            P_bonds = len(gene) - 1
            tot_n += Nn
            gene_meanlength = get_mean_length(tot_n, len(self.gene_list))

            # Gene and genome weight
            gene_weight = get_weight(DNA_cur_dict) - (P_bonds * self.Metabolites['H2O']['MW'])
            genome_weight += gene_weight # Weight of the genome (Da)


            '''---------------------'''
            '''Gene Formation energy'''
            dGf_gene1 = self.get_dGf_na(gene, DNA_cur_dict,'DNA','block')

            '''---------------------'''
            '''Gene Reaction energy'''
            dGr_gene1 = self.get_dGr_na(gene, DNA_cur_dict, dGf_gene1, 'minimum')

            DNA_dGr_list1.append(dGr_gene1)

        '''Molarity'''
        nucleotide_av_mw = genome_weight / tot_n # Average nucleotide weight (Da)
        av_gene_w = genome_weight / len(self.gene_list) # Average gene weight (Da)
        DNA_moles = self.cell.get_DNA_weight() / nucleotide_av_mw # Moles of DNA
        DNA_molarity = DNA_moles / self.cell.get_volume()



        '''Reqction quotient'''
        R = 0.008314472 #Gas constant KJ mol K

        for gene1, gene_dGr1 in zip(self.gene_list, DNA_dGr_list1):
            gene_lnQ1 = get_lnQ(DNA_molarity, gene, len(self.gene_list), self.Nucleic_Acids['DNA'], self.cell.get_model())

            '''Molar Gibbs Energy of each gene'''
            DNA_dG1 = (gene_dGr1 + (R * self.T * gene_lnQ1))
            DNA_reactionGs1.append(DNA_dG1)

        '''Genome energy (KJ/g)'''
        genome_energy1 = get_omic_energy(DNA_reactionGs1, av_gene_w)

        return genome_energy1


    def genome_synthesis2(self):
        """
        Calculate and return the free energy of synthesising one dry g of the
        genome of this instance's Cell.
        """

        #Nucleotide count for each gene
        #--- Number of nucleotides on each strand and total weight---
        '''-----DNA-----'''
        tot_n = 0
        genome_weight = 0
        #---Energy of each strand (Chain method )--
        DNA_dGf_list = []
        DNA_dGr_list2 = []
        DNA_reactionGs2 = []

        for gene in self.gene_list:

            '''General characteristics of the gene'''
            #Nucleotide count and gene lenght
            DNA_cur_dict = deepcopy(self.Nucleic_Acids['DNA'])
            DNA_n_count = Counter(gene)

            for n, dictionary in DNA_cur_dict.items():
                dictionary['count'] = DNA_n_count.get(n ,0)

            Nn = len(gene)
            P_bonds = len(gene) - 1
            tot_n += Nn
            gene_meanlength = get_mean_length(tot_n, len(self.gene_list))

            # Gene and genome weight
            gene_weight = get_weight(DNA_cur_dict) - (P_bonds * self.Metabolites['H2O']['MW'])
            genome_weight += gene_weight # Weight of the genome (Da)


            '''---------------------'''
            '''Gene Formation energy'''
            dGf_gene2 = self.get_dGf_na(gene, DNA_cur_dict,'DNA','block2')

            '''---------------------'''
            '''Gene Reaction energy'''
            dGr_gene2 = self.get_dGr_na(gene, DNA_cur_dict, dGf_gene2, 'minimum')

            DNA_dGr_list2.append(dGr_gene2)

        '''Molarity'''
        nucleotide_av_mw = genome_weight / tot_n # Average nucleotide weight (Da)
        av_gene_w = genome_weight / len(self.gene_list) # Average gene weight (Da)
        DNA_moles = self.cell.get_DNA_weight() / nucleotide_av_mw # Moles of DNA
        DNA_molarity = DNA_moles / self.cell.get_volume()

        '''Reqction quotient'''
        R = 0.008314472 #Gas constant KJ mol K
        for gene2, gene_dGr2 in zip(self.gene_list, DNA_dGr_list2):
            gene_lnQ2 = get_lnQ(DNA_molarity, gene, len(self.gene_list), self.Nucleic_Acids['DNA'], self.cell.get_model())

            '''Molar Gibbs Energy of each gene'''
            DNA_dG2 = (gene_dGr2 + (R * self.T * gene_lnQ2))
            DNA_reactionGs2.append(DNA_dG2)

        '''Genome energy (KJ/g)'''
        genome_energy2 = get_omic_energy(DNA_reactionGs2, av_gene_w)

        return genome_energy2

    def transcriptome_synthesis(self):
        """
        Calculate and return the free energy of synthesising one dry g of the
        transcriptome of this instance's Cell.
        """

        RNA_tot_n = 0
        transcriptome_weight = 0
        #---Energy of each strand (Chain method )--
        RNA_dGf_list = []
        RNA_dGr_list = []
        RNA_reactionGs = []
        for m in self.mRNA_list:

            '''General characteristics of the transcript'''
            #Nucleotide count and Trasncript lenght

            RNA_cur_dict = deepcopy(self.Nucleic_Acids['RNA'])
            RNA_n_count = Counter(m)

            for RNA_n, RNA_dictionary in RNA_cur_dict.items():
                RNA_dictionary['count'] = RNA_n_count.get(RNA_n ,0)

            RNA_Nn = len(m)
            RNA_P_bonds = len(m) - 1
            RNA_tot_n += RNA_Nn
            RNA_meanlength = get_mean_length(RNA_tot_n, len(self.mRNA_list))

            #Transcript and nucleotides Weight
            RNA_weight = get_weight(RNA_cur_dict) - (RNA_P_bonds * self.Metabolites['H2O']['MW'])
            transcriptome_weight += RNA_weight # Weight of the transcriptome (Da)


            '''---------------------'''
            '''Transcript Formation energy'''
            dGf_transcript = self.get_dGf_na(m, RNA_cur_dict,'RNA','chain')
            RNA_dGf_list.append(dGf_transcript)

            '''---------------------'''
            '''Transcript Reaction energy'''
            dGr_transcript = self.get_dGr_na(m, RNA_cur_dict, dGf_transcript, 'minimum')
            RNA_dGr_list.append(dGr_transcript)

        '''Molarity'''
        RNA_nucleotide_av_mw = transcriptome_weight / RNA_tot_n #Nucleotide average weight (Da)
        av_transcript_w = transcriptome_weight / len(self.mRNA_list) #Average transcript weight (Da)
        RNA_moles = self.cell.get_RNA_weight() / RNA_nucleotide_av_mw #mol
        RNA_molarity = RNA_moles / self.cell.get_volume()

        '''Reqction quotient'''
        R = 0.008314472 #Gas constant KJ mol K
        for m, RNA_dGr in zip(self.mRNA_list, RNA_dGr_list):
            m_lnQ = get_lnQ(RNA_molarity, m, len(self.mRNA_list), self.Nucleic_Acids['RNA'], self.cell.get_model())


            '''Molar Gibbs Energy of each transcript'''
            RNA_dG = (RNA_dGr + (R * self.T * m_lnQ))
            RNA_reactionGs.append(RNA_dG)

        '''Transcriptome energy (KJ/g)'''
        transcriptome_energy = get_omic_energy(RNA_reactionGs, av_transcript_w)

        # print('Transcriptome energy', transcriptome_energy, 'KJ/g at: ', T, 'K')

        return transcriptome_energy
