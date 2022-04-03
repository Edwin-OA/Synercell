"""
Cell, a class for storing all relevant information about a biological cell
for use throughout Syncell.

Cells can be initialised manually or by reading from the CellData csv file.
To add new organisms, try out adding rows to that file and using the relevant
keywords!

"""
from pandas import DataFrame, read_csv
import pandas as pd
import sys, os

from time import sleep


class Cell:

    ''' Constructors '''

    def __init__(self, model, cell_wet_weight, cell_dry_weight, DNA_weight, RNA_weight, proteins_weight, lipids_weight, glycans_weight, ions_weight, volume, Genome, Proteome, Genome_energy):
        self.model = model
        self.cell_wet_weight = cell_wet_weight
        self.cell_dry_weight = cell_dry_weight
        self.DNA_weight = DNA_weight
        self.RNA_weight = RNA_weight
        self.proteins_weight = proteins_weight
        self.lipids_weight = lipids_weight
        self.glycans_weight = glycans_weight
        self.ions_weight = ions_weight
        self.volume = volume
        self.Genome = Genome
        self.Proteome = Proteome
        self.Genome_energy = Genome_energy


    @classmethod
    def from_data(cls, celltype, Genome, Proteome, Genome_energy = 0.):
        """
        Construct a cell from tabulated data.

        Current celltypes include 'EC', 'SC', 'M'

        :param: celltype - keyword for the desired species.
        :param: Genome - filename of a genomic fasta file.
        :param: Proteome - filename of a proteomic fasta file.
        :param: Genome_energy - energy of genome (default 0. before calculation)
        """

        df = pd.read_csv(os.path.dirname(__file__)+'/CellData.csv')
        _df = df.loc[df['Keyword'] == celltype]
        for index, row in _df.iterrows():
            return cls(row['Keyword'],
              row['cell_wet_grams [g]'],
              row['cell_dry_grams [g]'],
              row['DNA_dry_grams [g]'],
              row['RNA_dry_grams [g]'],
              row['Proteins_dry_grams [g]'],
              row['Lipids_dry_grams [g]'],
              row['Glycans_dry_grams [g]'],
              row['ions_dry_grams [g]'],
              row['cell_intracellular_litres [L]'],
              Genome,
              Proteome,
              Genome_energy)


    """ get methods """

    def get_model(self):
        return self.model

    def get_cell_wet_weight(self):
        return self.cell_wet_weight

    def get_cell_dry_weight(self):
        return self.cell_dry_weight

    def get_DNA_weight(self):
        return self.DNA_weight

    def get_RNA_weight(self):
        return self.RNA_weight

    def get_proteins_weight(self):
        return self.proteins_weight

    def get_lipids_weight(self):
        return self.lipids_weight

    def get_glycans_weight(self):
        return self.glycans_weight

    def get_ions_weight(self):
        return self.ions_weight

    def get_volume(self):
        return self.volume

    def get_Genome(self):
        return self.Genome

    def get_Proteome(self):
        return self.Proteome

    def get_Genome_energy(self):
        return self.Genome_energy
