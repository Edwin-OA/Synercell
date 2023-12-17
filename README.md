# Synercell

Synercell is a computational tool designed to estimate the minimum energy required for synthesising a cell from its components. Targeted towards bioinformaticians, biologists, and researchers in synthetic biology, Synercell leverages advanced algorithms to provide insights into the bioenergetics of cells.


<h2>Installation and Setup </h2>

Prerequisites:
Python 3.X
Biopython, Matplotlib, ProgressBar, tqdm



<h2>Step-by-Step Installation: </h2>

1. Install Miniconda: Download from [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and follow the installation instructions.

2. Set Up Conda Environment:
```
conda config --append channels conda-forge
conda create -n rkt reaktoro
conda activate rkt

```

3. Install Dependencies:
```
pip install biopython matplotlib progressbar2 tqdm
```


<h2> Usage Instructions </h2>
Run Syncell using the command:

```
python _syncell.py.
```

The tool will prompt you for the following inputs:

**Cell Type:** Choose from bacterial, yeast, or mammalian.<br />
**Genome Sequence:** Provide the path to the genome sequence file.\* <br />
**Protein Sequence:** Provide the path to the protein sequence file.\* <br />
**Temperature:** Specify the temperature for the calculation. ± <br />

\* It is recommended to have the sequence file in  in the same directory as the synercell.py module. FASTA or FNA files are preferred. <br />
± The temperature range supported is 275 K  to 400 K.

<h2> Workflow: </h2>
To run the program: 

```
python _syncell.py 
```

Follow on-screen prompts to input, cell type, sequences and temperature.

**-----Welcome to Syncell-----**
 
This program will help you obtain the energetic cost of synthesising a cell by calculating the energy of the genome, transcriptome, proteome and membrane:
 
 **-Model size selection-**
To start, choose between one of the three model cell sizes available: E. coli [EC], S. cerevisiae [SC] or average mammalian [M]

> Choose one of the available cell sizes available, in this case, bacteria:

EC / SC / M: **EC**
 
**-Temperature selection-**
Do you want to obtain the information of a single temperature point from the range 275 K - 400 K?
[yes/no]:  **yes**

> If you select yes, the program will ask for one temperature to display the results. In either case, it will display a graph in the end with all the temperature points from 275 K to 400 K

Type the temperature (in kelvin) you would like to print results for: **298**
 
 **-Sequence input-**
You can choose between calculating the Proteome energy, the Genome, trancriptome and membrane energy or all.

> You can choose between calculating the energy of the genome, transcriptome and lipid bilayer (N), the proteome and lipid bilayer (P), or the genome, transcriptome, proteome and lipid bilayer (A)

Nucleic Acids [N] / Proteins [P] / All [A]: **A**

> After choosing the type of calculation, the program will ask for the genome sequence on a fasta file, the proteome sequence on a fasta file or both.

First, enter the name of a DNA fasta file (e.g. Ecoli_genome.fasta)
DNA file name: **ecoli_genome.fasta**

Now, enter the name of a fasta file containing one or more protein sequences (e.g. Ecoli_proteome.fasta)

Protein file name: **ecoli_proteome.fasta**


--genome genome.fasta --protein proteins.fasta


Interpreting Results: The output will be displayed in the console in specified formats.


<h2> Output and Interpretation </h2>

Format: Results are outputted as five plots, containing the energy necessary to synthesise the cell's components in Joules. The plots will appear in the following order:

1. Genome and transcriptome
2. Proteome
3. Lipid bilayer
4. All cell components comparison
5. Cell 




<h2> Methodology </h2>

Synercell applies the Group Contribution Algorithm (GCA) and Gibbs Free Energy of biomolecules' building blocks to estimate the standard Gibbs free energy of formation for proteins, DNA, RNA, and lipids. For a detailed scientific explanation, refer to our paper.


<h2> Contributing </h2>

We welcome contributions! Please refer to our contribution guidelines for how to contribute or get in touch: edwin.oa@icloud.com.

