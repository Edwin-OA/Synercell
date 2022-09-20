# Syncell
A program to create a virtual cell and calculate its energy of synthesis.


Getting ready to use Syncell

To run Syngen and Synpro contained in Syncell is necessary to install the reaktoro environment and the biopython, matplotlib, progressbar and tqdm packages in case you don't have them already. To do this simply follow the next steps:

Reaktoro

1. In the terminal choose the work carpet where the programs and your files are contained
2. Instal miniconda from the following link https://docs.conda.io/en/latest/miniconda.html
3. Download the 64-bit bash file
4. In the terminal window type 'bash' and drag the downloaded file next to bash in the terminal window
5. Follow the instructions to install miniconda and accept the terms and conditions
6. Next, add the conda-forge channels by typing in the terminal: conda config --append channels conda-forge
7. Once done this, install reaktoro by typing: conda install reaktoro=1
8. To create the environment for reaktoro type: conda create -n rkt reaktoro
9. Before running Syncell is necessary to activate reaktoro by typing: conda activate rkt or conda activate (rkt path) Biopython
10. In the terminal window type: pip install biopython Matplotlib
11. In the terminal type: python -m pip install -U matplotlib ProgressBar
12. Similarly, type: pip install progressbar2 ProgressBar
13 Finally, type: pip install tqdm
