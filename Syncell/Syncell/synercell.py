"""Synercell. Version E.2 - 05.01.24"""

import sys, os, math, statistics, re
sys.path.append('Syncell/')
import os.path
import pandas as pd
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

# --- Welcome message and instructions ---

def main():

    Syncell_list = []
    TK = np.linspace(275,400, num = 3)

    print('\n----- Welcome to Synercell -----\n')
    print('Synercell is designed to assist you in calculating the energy required to synthesise a cell at different temperatures. Here is how to get started:\n')
    print('\n1. Select a Cell Size\n')
    print('\nTo use a concentration pool, begin by choosing a cell size category from one of the three available options:\n')
    print('Bacteria: Represented by Escherichia coli Code: EC')
    print('Yeast: Represented by Saccharomyces cerevisiae Code: SC')
    print('Mammalian: An average mammalian cell Code: M')
    print(' ')

    concentration_pool = get_concentration_pool()

    print(' ')
    print('\n2. Enter a Cell Code\n') 
    print(' ')
    print('Input the code corresponding to the type of cell for which you want to estimate the energy cost. You can use one of the standard codes (EC, SC, M) or')
    print('You may also enter a custom code from the CellData.csv file.\n')

    cell_code = get_cell_code()

    print(' ')
    print('\n3. Temperature selection')
    print('If you wish to obtain the information of a single temperature point from the range 275 K - 400 K, select the temperature you are interested in:')

    selected_temp, TK = temperature_selection(TK)

    print(' ')
    print('\n4. Tool selection')
    print('You can choose between calculating the Proteome energy, the Genome, trancriptome and lipid bilayer energy or al. lThe programm will always calculate the lipid bilayer energey by default')

    tool, gname, pname = tool_input()

    # Create a Cell object using the chosen parameters
    cell = Cell.from_data(cell_code, gname, pname, tool, Genome_energy=0.)
   
    # Define temperature range based on the selected temperature
    TK = define_temperature_range(selected_temp)

    # Calculate the energetic costs for the chosen parameters
    calculate_energetic_costs(tool, cell, TK)

    # Additional processing or output can be added here as needed.

# --- Function definitions ---

def get_concentration_pool():
    while True:
        concentration_pool = input('Enter cell size (EC / SC / M): ').lower()
        if concentration_pool in ['ec', 'bacteria']:
            return 'EC', np.linspace(275,400, num=5)
        elif concentration_pool in ['sc', 'yeast']:
            return 'SC', np.linspace(275,400, num=5)
        elif concentration_pool in ['m', 'mammal']:
            return 'M', np.linspace(275,400, num=2)
        elif concentration_pool == 'other':
            cell_code = input('Enter two-character cell code (e.g., 3A): ')
            return cell_code, np.linspace(275,400, num=3)
        elif concentration_pool == 'q':
            quit()

def get_cell_code():
    while True:
        cell_code = input('Enter cell code (EC / SC / M / Other): ').lower()
        if cell_code in ['ec', 'bacteria']:
            return 'EC'
        elif cell_code in ['sc', 'yeast']:
            return 'SC'
        elif cell_code in ['m', 'mammal']:
            return 'M'
        elif cell_code == 'other':
            return input('Enter two-character cell code (e.g., 3A): ')
        elif cell_code == 'q':
            quit()


def temperature_selection(TK):
    while True:
        temp = input('Would you like to specify a temperature? [yes/no]: ')
        if temp.lower() == 'no':
            return None, TK
        if temp.lower() == 'yes':
            while True:
                try:
                    stemperature = int(input("Enter the temperature (in Kelvin) for analysis: "))
                    if 275 <= stemperature <= 400:
                        TK = np.append(TK, stemperature)
                        TK = np.sort(TK)
                        selected_temp = np.nonzero(TK == stemperature)[0][0]
                        return selected_temp, TK
                    else:
                        print("Temperature out of range, try choosing between 275 and 400")
                except ValueError:
                    print("Please enter a valid number. Alternatively, type 'q' to exit.")
        if temp.lower() == 'q':
            quit()

def tool_input():
    while True:
        input_tool = input('Nucleic Acids [N] / Proteins [P] / All [A]: ' )
        if input_tool == "N" or input_tool == 'n' or input_tool == 'Nucleic Acids' or input_tool == 'nucleic acids':
            tool = 'Syngen'
            print('Enter the name of a DNA fasta file (e.g. Ecoli_genome.fasta)')
            gname = input('DNA file name: ')
            if os.path.exists(gname) == False:
                print('The file is not in the current directory')
                continue
            pname = []
            break
        if input_tool == "P" or input_tool == 'p' or input_tool == 'Proteins' or input_tool == 'proteins':
            tool = 'Synpro'
            print('Enter the name of a fasta file containing one or more protein sequences (e.g. Ecoli_proteome.fasta)')
            pname = input('Protein file name: ')
            gname = []
            if os.path.exists(pname) == False:
                print('The file is not in the current directory')
                continue
            break
        if input_tool == 'All' or input_tool == 'all' or input_tool == 'A' or input_tool == 'a':
            tool = 'both'
            print('Now, enter the name of a DNA fasta file (e.g. Ecoli_genome.fasta)')
            gname = input('DNA file name: ')
            if os.path.exists(gname) == False:
                print('The file is not in the current directory')
                continue
            print('Now, enter the name of a fasta file containing one or more protein sequences (e.g. Ecoli_proteome.fasta)')
            while True:
                pname = input('Protein file name: ')
                if pname == 'q':
                    quit()
                if os.path.exists(pname) == False:
                    print('The file is not in the current directory')
                    continue
                if os.path.exists(pname) == True:
                    break
            break
        if input_tool == 'q':
            quit()


# def define_temperature_range(selected_temp):
#     # Logic to define the temperature range based on the selected temperature
#     pass

# def calculate_energetic_costs(tool, cell, TK):
#     # Logic to calculate energetic costs
#     pass
    TK = np.linspace(275,400, num=24)

    # if temp == 'yes':
    #     TK=np.append(TK,stemperature)
    #     TK=np.sort(TK)
    #     selected_temp= np.nonzero(TK == stemperature)[0][0]


    '''----Syngen implementation---'''

    '''Genome'''
    dG_genome_T = []
    '''Transcriptome'''
    dG_transcriptome_T = []
    '''Proteome'''
    dG_proteome_T = []
    bio_dG_proteome_T = []
    '''Membrane'''
    dG_membrane_T = []
    '''Cell's energy'''
    Cell_energy_t = []

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

            '''Transcriptome calculations'''
            transcriptome_energy = SG.transcriptome_synthesis()

            '''Genome list'''
            dG_genome_T.append(genome_energy)

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

        ax.plot(TK, dG_transcriptome_T, label='Transcriptome Chain Energy', color='yellow', marker='1', linewidth=2)

        ax.set_ylabel(r'Energetic cost [KJ per gram]', fontsize=14)
        ax.set_xlabel('Temperature [K]', fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14)

        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.show()


    '''Plotting proteome energy'''
    if tool == 'Synpro' or tool == 'both':

        fig = plt.figure(figsize = (7,5))
        ax= fig.add_subplot(111)
        plt.title('Molar Gibbs Energy calculated at two conditions')
        ax.plot(TK, dG_proteome_T, label='Proteome energy using standard energies', color='lawngreen', marker='v', linewidth=3)
        # ax.plot(TK, bio_dG_proteome_T, label='Proteome energy using biological standard energies', color='steelblue', marker='^', linewidth=2)

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
        G1 = ax.plot(TK, dG_genome_T, label='Genome Chain Energy', color='crimson', marker='1', linewidth=2)

        '''Transcriptome'''
        T1 = ax.plot(TK, dG_transcriptome_T, label='Transcriptome Chain Energy', color='yellow', marker='1', linewidth=2)

        '''Proteome'''
        P1 = ax.plot(TK, dG_proteome_T, label='Proteome energy using standard energies', color='lawngreen', marker='v', linewidth=3)


        '''Membrane'''
        M1 = ax.plot(TK, dG_membrane_T, label='Membrane energy', color='purple', marker='v', linewidth=3)

        ax.set_ylabel(r'Energetic cost [kJ per dry g]', fontsize=14)
        ax.set_xlabel('Temperature [K]', fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=14)
        plt.legend(fontsize=14)
        plt.tight_layout()
        plt.show()

        results = {'Organism': cell_model,'Temperature': TK, 'Genome_chain':dG_genome_T, 'Genome_B1':dG_genome_T1, 'Genome_B2':dG_genome_T2, 'RNA_chain':dG_transcriptome_T,  'thermo_prot': dG_proteome_T, 'bio_prot': bio_dG_proteome_T, 'membrane': dG_membrane_T}
        df = pd.DataFrame(results)
        df.to_csv('Results.csv')

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
                print('Genome energy by block', "%.3f" % genome_energy_t1, 'KJ/g')
                print('Genome energy by block 2', "%.3f" % genome_energy_t2, 'KJ/g')

                print('Transcriptome energy', "%.3f" % transcriptome_energy_t, 'KJ/g')
            if tool == 'Synpro' or tool == 'both':
                print('Proteome energy', "%.3f" % proteome_energy_t, 'KJ/g')
            if tool == 'Synpro' or tool == 'Syngen' or tool == 'both':
                print('Membrane energy', "%.3f" % membrane_energy_t, 'KJ/g')
                Syncell_t = (sum(Syncell_list) * cell.cell_dry_weight) * 1000
                print('Energy to synthesie one single cell', Syncell_t, 'J per cell at', stemperature, 'K')

# --- Main execution ---

if __name__ == '__main__':
    main()