"""Syncell. Version C.1 - 16.10.21"""

import sys, os, math, statistics, re
sys.path.append('Syncell/')
import os.path
from multiprocessing import Process,Queue,Pipe
import subprocess
import numpy as np
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.Seq import transcribe
import matplotlib.pyplot as plt
from copy import deepcopy
from collections import Counter
# from BioMolecule import BioMolecule


# time and progressbars
from time import sleep
import time
import progressbar
from tqdm import tqdm, trange
from progressbar import AnimatedMarker, Bar, BouncingBar, ETA, \
    AdaptiveETA, FileTransferSpeed, FormatLabel, Percentage, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer, UnknownLength


from Syncell import Syngen, Synpro, Synpho, Cell


# --- functions ---

def main():

    Syncell_list = []

    '''----Input block---'''

    print(' ')
    print('-----Welcome to Syncell-----')
    print(' ')
    print('This program will help you obtain the energetic cost of synthesising a cell by calculating the energy of the genome, transcriptome, proteome and membrane:')
    print(' ')
    print(' -Model size selection-')
    print('To start, choose between one of the three model cell sizes available: E. coli [EC], S. cerevisiae [SC] or average mammalian [M]')
    while True:
        input_cell_model= input('EC / SC / M: ' )
        if input_cell_model == 'EC' or input_cell_model == 'ec' or input_cell_model == 'bacteria':
            cell_model = 'EC'
            TK = np.linspace(275,400, num=5)
            break
        if input_cell_model == 'SC' or input_cell_model == 'sc' or input_cell_model == 'yeast':
            cell_model = 'SC'
            TK = np.linspace(275,400, num=5)
            break
        if input_cell_model == 'M' or input_cell_model == 'm' or input_cell_model == 'mammal':
            cell_model = 'M'
            TK = np.linspace(275,400, num=2)
            break
        if input_cell_model == 'q':
            quit()
        else:
            print("We could not get that. Try typing only the key words in capital letters.")
            print("Alternatively type 'q' to exit")
    print(' ')
    print(' -Temperature selection-')
    print('Do you want to obtain the information of a single temperature point from the range 275 K - 400 K?')
    while True:
        temp = input('[yes/no]: ')
        if temp == 'no':
            stemperature = 0
            selected_temp =[]
            break
        if temp == 'yes':
            while True:
                try:
                    stemperature = int(input("Type the temperature (in kelvin) you would like to print results for:"))
                except ValueError:
                    print("Try using only numbers, alternatevely type 'q' to exit the program")
                    continue
                if stemperature < 249 or stemperature > 400:
                    print("Temperature out of range, try choosing between 275 and 400")
                    continue
                if stemperature > 250 or stemperature < 400:
                    TK = np.append(TK, stemperature)
                    TK = np.sort(TK)
                    selected_temp= np.nonzero(TK == stemperature)[0][0]
                    break
                if stemperature == 'q':
                    quit()
            break
        if temp == 'q':
            quit()
            break

    print(' ')
    print(' -Sequence input-')
    print('You can choose between calculating the Proteome energy, the Genome, trancriptome and membrane energy or all')
    while True:
        input_tool = input('Nucleic Acids [N] / Proteins [P] / All [A]: ' )
        if input_tool == "N" or input_tool == 'n' or input_tool == 'Nucleic Acids' or input_tool == 'nucleic acids':
            tool = 'Syngen'
            print('Enter the name of a DNA fasta file (e.g. Ecoli_genome.fasta)')
            gname = input('DNA file name: ')
            if os.path.exists(gname) == False:
                print('The file is not in current directory')
                continue
            pname = []
            break
        if input_tool == "P" or input_tool == 'p' or input_tool == 'Proteins' or input_tool == 'proteins':
            tool = 'Synpro'
            print('Enter the name of a fasta file containing one or more protein sequences (e.g. Ecoli_proteome.fasta)')
            pname = input('Protein file name: ')
            gname = []
            if os.path.exists(pname) == False:
                print('The file is not in current directory')
                continue
            break
        if input_tool == 'All' or input_tool == 'all' or input_tool == 'A' or input_tool == 'a':
            tool = 'both'
            print('First, enter the name of a DNA fasta file (e.g. Ecoli_genome.fasta)')
            gname = input('DNA file name: ')
            if os.path.exists(gname) == False:
                print('The file is not in current directory')
                continue
            print('Now, enter the name of a fasta file containing one or more protein sequences (e.g. Ecoli_proteome.fasta)')
            while True:
                pname = input('Protein file name: ')
                if pname == 'q':
                    quit()
                if os.path.exists(pname) == False:
                    print('The file is not in current directory')
                    continue
                if os.path.exists(pname) == True:
                    break
            break
        if input_tool == 'q':
            quit()

    '''----Cell definition---'''

    cell = Cell.from_data(cell_model, gname, pname, Genome_energy=0.)

    '''----Syngen implementation---'''

    '''Genome'''
    dG_genome_T = []
    dG_genome_T1 = []
    dG_genome_T2 = []
    '''Transcriptome'''
    dG_transcriptome_T = []
    '''Proteome'''
    dG_proteome_T = []
    bio_dG_proteome_T = []
    '''Membrane'''
    dG_membrane_T = []
    '''Cell's energy'''
    Cell_energy_t = []

    TK = np.linspace(275,400, num=24)
    if temp == 'yes':
        TK=np.append(TK,stemperature)
        TK=np.sort(TK)
        selected_temp= np.nonzero(TK == stemperature)[0][0]

    '''Nucleic Acids'''
    if tool == 'Syngen' or tool == 'both':
        # Generate a Syngen object
        SG = Syngen(TK[0], cell)
        widget1 = ['Calculating Gibbs Energy of the nucleic acids', progressbar.Bar('§'), '(', progressbar.Percentage(), ' complete)',
        ' [', progressbar.Timer(), '] '
        ' (', progressbar.ETA(), ') ',
        ]

        firststep = True
        for T in progressbar.progressbar(TK, widgets = widget1):

            if not firststep:
                SG.update_T(T)
            else:
                firststep = False

            '''Genome calculations'''
            genome_energy = SG.genome_synthesis()
            genome_energy1 = SG.genome_synthesis1()
            genome_energy2 = SG.genome_synthesis2()

            '''Transcriptome calculations'''
            transcriptome_energy = SG.transcriptome_synthesis()

            '''Genome list'''
            dG_genome_T.append(genome_energy)
            dG_genome_T1.append(genome_energy1)
            dG_genome_T2.append(genome_energy2)

            '''Transcriptome list'''
            dG_transcriptome_T.append(transcriptome_energy)

        if temp == 'yes':
            genome_energy_t = dG_genome_T[selected_temp]
            transcriptome_energy_t = dG_transcriptome_T[selected_temp]
            Syncell_list.extend([genome_energy_t, transcriptome_energy_t])


    '''Proteins'''
    if tool == 'Synpro' or tool == 'both':
        # Generate a Synpro object
        SP = Synpro(TK[0], cell)
        widget2 = ['Calculating Gibbs Energy of the proteins', progressbar.Bar('|'), '(', progressbar.Percentage(), ' complete)',
            ' [', progressbar.Timer(), '] '
            ' (', progressbar.ETA(), ') ']


        firststep = True
        for T in progressbar.progressbar(TK, widgets=widget2):

            if not firststep:
                SP.update_T(T)
            else:
                firststep = False

            '''Proteome calculations'''
            proteome_energy, proteome_bio_energy = SP.proteome_synthesis()

            '''Proteome list'''
            dG_proteome_T.append(proteome_energy)
            bio_dG_proteome_T.append(proteome_bio_energy)

        if temp == 'yes':
            proteome_energy_t = bio_dG_proteome_T[selected_temp]
            Syncell_list.append(proteome_energy_t)

    '''Phospholipids'''
    if tool == 'Syngen' or tool == 'Synpro' or tool == 'both':
        # Generate a Synpho object
        SPh = Synpho(TK[0], cell)
        widget3 = ['Calculating Gibbs Energy of the membrane', progressbar.Bar('|'), '(', progressbar.Percentage(), ' complete)',
        ' [', progressbar.Timer(), '] '
        ' (', progressbar.ETA(), ') ']

        firststep = True
        for T in progressbar.progressbar(TK, widgets=widget3):

            if not firststep:
                SPh.update_T(T)
            else:
                firststep = False

            '''Membrane calculations'''
            membrane_energy = SPh.membrane_synthesis()

            '''Membrane list'''
            dG_membrane_T.append(membrane_energy)

        if temp == 'yes':
            membrane_energy_t = dG_membrane_T[selected_temp]
            Syncell_list.append(membrane_energy_t)

    '''Cell'''
    if tool == 'both':
        Cell_e = [genome + transcriptome + proteome + membrane for genome, transcriptome, proteome, membrane in zip(dG_genome_T, dG_transcriptome_T, dG_proteome_T, dG_membrane_T)]

        for system_t in Cell_e:
            joules = (system_t * cell.cell_dry_weight) * 1000
            Cell_energy_t.append(joules)

    '''Plotting Genome Energy'''
    if tool == 'Syngen' or tool == 'both':
        fig = plt.figure(figsize = (7,5))
        ax= fig.add_subplot(111)

        plt.title('Molar Gibbs Energy of the genome and transcriptome')
        ax.plot(TK, dG_genome_T, label='Genome Chain Energy', color='crimson', marker='1', linewidth=2)
        ax.plot(TK, dG_genome_T1, label='Genome Block 1 Energy ', color='blue', marker='1', linewidth=2)
        ax.plot(TK, dG_genome_T2, label='Genome Block 2 Energy', color='green', marker='1', linewidth=2)

        ax.plot(TK, dG_transcriptome_T, label='Transcriptome Chain Energy', color='yellow', marker='1', linewidth=2)

        ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
        ax.set_xlabel('Temperature [K]', fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14)

        #ax.set_xlim(270, 400)
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.show()


    '''Plotting proteome energy'''
    if tool == 'Synpro' or tool == 'both':

        fig = plt.figure(figsize = (7,5))
        ax= fig.add_subplot(111)
        plt.title('Molar Gibbs Energy calculated at two conditions')
        ax.plot(TK, dG_proteome_T, label='Proteome energy using standard energies', color='lawngreen', marker='v', linewidth=3)
        ax.plot(TK, bio_dG_proteome_T, label='Proteome energy using biological standard energies', color='steelblue', marker='^', linewidth=2)

        ax.set_ylabel(r'Energetic cost [kJ per dry g]', fontsize=14)
        ax.set_xlabel('Temperature [K]', fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_xlim(270, 400)
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.show()


    ''''Plotting membrane energy'''
    if tool == 'Syngen' or tool == 'Synpro' or tool == 'both':

        fig = plt.figure(figsize = (7,5))
        ax= fig.add_subplot(111)
        plt.title('Molar Gibbs Energy of the membrane')
        ax.plot(TK, dG_membrane_T, label='Membrane energy', color='purple', marker='v', linewidth=3)

        ax.set_ylabel(r'Energetic cost [kJ per dry g]', fontsize=14)
        ax.set_xlabel('Temperature [K]', fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_xlim(270, 400)
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.show()

    '''Plotting All'''
    if tool == 'both':
        fig = plt.figure(figsize = (7,5))
        ax= fig.add_subplot(111)
        plt.title('Molar Gibbs Energy of the cell structures')
        '''Genome'''
        ax.plot(TK, dG_genome_T, label='Genome Chain Energy', color='crimson', marker='1', linewidth=2)
        ax.plot(TK, dG_genome_T1, label='Genome Block 1 Energy ', color='blue', marker='1', linewidth=2)
        ax.plot(TK, dG_genome_T2, label='Genome Block 2 Energy', color='green', marker='1', linewidth=2)
        '''Transcriptome'''
        ax.plot(TK, dG_transcriptome_T, label='Transcriptome Chain Energy', color='yellow', marker='1', linewidth=2)
        '''Proteome'''
        ax.plot(TK, dG_proteome_T, label='Proteome energy using standard energies', color='lawngreen', marker='v', linewidth=3)
        ax.plot(TK, bio_dG_proteome_T, label='Proteome energy using biological standard energies', color='steelblue', marker='^', linewidth=2)
        '''Membrane'''
        ax.plot(TK, dG_membrane_T, label='Membrane energy', color='purple', marker='v', linewidth=3)

        ax.set_ylabel(r'Energetic cost [kJ per dry g]', fontsize=14)
        ax.set_xlabel('Temperature [K]', fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14)
        # ax.set_xlim(270, 400)
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.show()

    ''''Plotting cell energy'''
    if tool == 'both':
        fig = plt.figure(figsize = (7,5))
        ax= fig.add_subplot(111)
        plt.title('Molar Gibbs Energy of the cell')
        ax.plot(TK, Cell_energy_t, label='Cell energy', color='green', marker='v', linewidth=3)

        ax.set_ylabel(r'Energetic cost [Joules per Cell]', fontsize=14)
        ax.set_xlabel('Temperature [K]', fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_xlim(270, 400)
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.show()


        '''Final results'''
        if temp == 'yes':
            if tool == 'Syngen' or tool == 'both':
                print('')
                print('--------------------------')
                print('The energy required to synthesise the components of the cell at', stemperature , 'K is:')
                print('Genome energy', "%.3f" % genome_energy_t, 'KJ/g')
                print('Transcriptome energy', "%.3f" % transcriptome_energy_t, 'KJ/g')
            if tool == 'Synpro' or tool == 'both':
                print('Proteome energy', "%.3f" % proteome_energy_t, 'KJ/g')
            if tool == 'Synpro' or tool == 'Syngen' or tool == 'both':
                print('Membrane energy', "%.3f" % membrane_energy_t, 'KJ/g')
                #Equation 2 from Cell Biosynthesis Energetics
                Syncell_t = (sum(Syncell_list) * cell.cell_dry_weight) * 1000
                print('Energy to synthesie one single cell', Syncell_t, 'J per cell at', stemperature, 'K')


if __name__ == '__main__':
    main()
