# Synercell

This file provides guidance to humans and AI agents when working with code in this repository.  It describes the purpose of the software, how to run it, the required environment, and a high‑level overview of the architecture and data flow.  Use this document as a quick reference when automating tasks or exploring the code.

## Project Overview

Syncell is a computational biology tool for calculating the energetic cost of synthesising biological cells.  It computes Gibbs free energy for genome, transcriptome, proteome and membrane synthesis across temperature ranges (275–400 K) using thermodynamic data from the SUPCRT database via the Reaktoro library.  The tool supports multiple cellular models (E. coli, yeast, mammalian and synthetic minimal) and generates quantitative predictions of energy requirements under different conditions.

## Installation and setup

### Prerequisites

To run the software you will need Python 3.x and a handful of scientific Python libraries.  The core dependencies are Biopython for reading FASTA files, Matplotlib for plotting, and two progress‑bar libraries (`progressbar2` and `tqdm`) for providing feedback during long calculations.

### Step‑by‑step installation

#### Standard installation (Intel/x86_64 systems)

1. **Install Miniconda** – download the installer for your platform from the [Miniconda documentation](https://docs.conda.io/en/latest/miniconda.html) and follow the instructions to set it up.

2. **Configure the conda environment** – add the conda‑forge channel and create an environment named `rkt` containing Reaktoro version 1.  Then activate the environment:

   ```bash
   conda config --append channels conda-forge
   conda create -n rkt reaktoro=1 python=3.9
   conda activate rkt
   ```

3. **Install Python dependencies** – within the activated environment, install the required packages:

   ```bash
   pip install biopython matplotlib progressbar2 tqdm pandas
   ```

4. **Prepare input files** – place your genome and proteome FASTA files in the repository's root directory.  Genome files typically have a `.fasta` or `.fna` extension; proteome files usually end with `.faa` or `.fasta`.  The scripts expect to find these files in the current working directory.

#### Apple Silicon (M1/M2/M3) installation

Reaktoro v1 is not natively available for ARM64 (Apple Silicon) Macs.  However, it can be run using Rosetta 2 emulation:

1. **Install Miniconda** – follow the standard installation above.

2. **Accept conda Terms of Service** (if first time using conda):

   ```bash
   conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
   conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
   ```

3. **Create x86_64 environment** – use the `CONDA_SUBDIR` environment variable to force x86_64 architecture:

   ```bash
   conda config --append channels conda-forge
   CONDA_SUBDIR=osx-64 conda create -n rkt reaktoro=1 python=3.9 -y
   ```

4. **Configure environment for x86_64** – ensure the environment always uses x86_64:

   ```bash
   conda activate rkt
   conda config --env --set subdir osx-64
   ```

5. **Install Python dependencies**:

   ```bash
   pip install biopython matplotlib progressbar2 tqdm pandas
   ```

6. **Verify installation** – test that Reaktoro v1 works with the custom database:

   ```bash
   python test_reaktoro_v1.py
   ```

   You should see success messages indicating that the custom database loads correctly and all organic species are accessible.


## Commands

### Running the main programs

Two interactive scripts are provided.  The original command‐line interface can be started with:

```bash
python syncell.py
```

A newer version with an improved user interface is available via:

```bash
python synercell.py
```

Both programs are interactive and will prompt for:

1. **Cell model selection** – choose EC, SC, M or 3A for _E. coli_, _S. cerevisiae_, mammalian or JCVI‑syn3A models.
2. **Temperature selection** – optionally specify a single temperature in kelvin; otherwise a range from 275 K to 400 K will be used.
3. **Tool selection** – pick nucleic acids (N), proteins (P) or all (A) to compute genome/transcriptome, proteome or whole‑cell energies.
4. **Input FASTA files** – provide genome and/or proteome sequence files in FASTA format.

### Environment activation

Before running any Syncell programs, activate the conda environment:

```bash
conda activate rkt
```

**Note:** Reaktoro version 1 is required for compatibility with the custom SUPCRT database (`Syncell/supcrt07-organics_AABB.xml`).  See the **Installation and setup** section above for detailed setup instructions, including special instructions for Apple Silicon Macs.

If you haven't installed the environment yet, refer to the installation instructions at the top of this document.  For troubleshooting, see the **Troubleshooting** section at the end.

## Usage instructions

After activating the `rkt` environment and installing the dependencies, you can run one of the main programs.  The recommended entry point is `synercell.py`, which provides an improved user interface.  The older `syncell.py` script remains available for backwards compatibility.

```bash
python synercell.py   # recommended interface

python syncell.py     # legacy interface
```

When the program starts it will guide you through the following prompts:

1. **Cell model** – choose one of the predefined models: EC (bacterial), SC (yeast), M (mammalian) or 3A (JCVI‑syn3A).
2. **Temperature selection** – indicate whether you want results at a single temperature or across the full range.  If you answer “yes” to the single‑temperature option, you will be asked for a temperature in kelvin between 275 K and 400 K.  Otherwise the default range (275–400 K in 24 steps) will be used.
3. **Calculation mode** – specify which biomolecules to analyse: N for nucleic acids (genome, transcriptome and membrane), P for proteins (proteome and membrane), or A to compute energies for all components.
4. **Input files** – provide the names of the genome and/or proteome FASTA files when prompted.  If you select the N mode, only a genome file is required; P mode requires only a proteome file; A mode requires both.

Respond to each prompt with the appropriate value.  Once all inputs are provided, the program will calculate Gibbs energies for the selected components over the specified temperature(s) and display plots summarising the results.  When run in “All” mode, the numerical results for each temperature are also saved to a CSV file for further analysis.


## Architecture

### Core module structure

The `Syncell/` directory contains the main computational classes:

* **Cell** (`Cell.py`) – data class storing cell properties (weights, volumes and references to genome/proteome files).  It loads cell specifications from `CellData.csv` using keywords EC, SC, M and 3A.  Use `Cell.from_data(celltype, genome_file, proteome_file)` to construct a cell model from tabulated data.

* **Syngen** (`Syngen.py`) – calculates energetic costs for genome and transcriptome synthesis.  It builds temperature‑dependent dictionaries of nucleic acid thermodynamic properties.  Key methods include `genome_synthesis()`, `genome_synthesis1()`, `genome_synthesis2()` and `transcriptome_synthesis()`.  Use `update_T(T)` to refresh internal properties when the temperature changes.

* **Synpro** (`Synpro.py`) – computes energetic costs for proteome synthesis.  It manages amino‑acid thermodynamic dictionaries (20 standard amino acids plus backbone species).  The main method `proteome_synthesis()` returns energies in both chemical and biological standard states.  Group contribution functions are implemented via `get_dGf_P()` and `get_dGr_P()`.

* **Synpho** (`Synpho.py`) – estimates energetic costs for membrane (phospholipid bilayer) synthesis using a group contribution approach for complex lipids such as POPC.  Its `membrane_synthesis()` method builds lipid energies from metabolic precursors (glucose, pyruvate, serine).

* **BioMolecule** (`BioMolecule.py`) – wraps Reaktoro thermodynamic calculations, converting chemical standard Gibbs energies to biological standard (pH 7) using Alberty’s 1998 transformation equations.  It reads data from `supcrt07-organics_AABB.xml`.

* **GCAtools** (`GCAtools.py`) – provides utility functions for group contribution analysis.  `get_lnQ()` computes reaction quotients for non‑standard conditions, `get_omic_energy()` converts molar energies to per‑gram values, and `get_Metabolites_dict(T)` builds temperature‑dependent metabolite dictionaries.

### Data flow

The typical workflow within the scripts is:

1. The user selects a cell model, and the program loads parameters from `CellData.csv` via the `Cell` class.
2. Genome and proteome sequences are read via BioPython and stored in the `Cell` object.
3. A temperature range is defined (by default 275–400 K with 24 points, or a single temperature if specified).
4. For each temperature:
   - Instances of `Syngen`, `Synpro` and `Synpho` are created or updated for that temperature.
   - `BioMolecule` queries Reaktoro for thermodynamic data at the given temperature.
   - Each synthesis class calculates Gibbs energies per dry gram of biomolecule.
5. The results are plotted with matplotlib and can be saved to `Results.csv` when running in “All” mode.

### Key thermodynamic concepts

* **Chemical standard state:** Gibbs energies from the SUPCRT database represent the chemical standard state.
* **Biological standard state:** Energies are converted to pH 7 via the `biogibbs()` transformation in `BioMolecule`.
* **Group contribution:** Complex molecules (particularly lipids) are built from simple precursors via group contributions in `Synpho`.
* **Concentration correction:** Cell‑type‑specific intracellular concentrations (EC, SC, M) influence reaction quotients through `get_lnQ()`.

### Cell models in `CellData.csv`

`CellData.csv` defines four model organisms with measured dry weights for DNA, RNA, proteins, lipids, glycans, ions and intracellular volume:

| Keyword | Organism                          | Description                      |
|--------:|------------------------------------|----------------------------------|
| EC      | _Escherichia coli_                | Bacterial model, smallest cell   |
| SC      | _Saccharomyces cerevisiae_        | Yeast model, intermediate size   |
| M       | Mammalian (average)               | Average mammalian cell, largest  |
| 3A      | JCVI‑syn3A                        | Minimal synthetic genome organism|

## Input file requirements

### Genome files

Genome files must be in FASTA format containing DNA sequences.  They may consist of multiple chromosomes or contigs.  Example files include `Ecoli_genome.fasta`, `Scerevisiae_genome.fna` and `3A_genome.fna`.

### Proteome files

Proteome files should also be in FASTA format with single‑letter amino acid codes, one entry per protein.  Example files include `Ecoli_proteome.fasta`, `Scerevisiae_proteome.fasta` and `3A_proteome.faa`.  All input files should reside in the repository’s root directory.

## Output

Both scripts generate matplotlib figures showing:

* Genome/transcriptome energy versus temperature
* Proteome energy (chemical standard vs. biological standard) versus temperature
* Membrane energy versus temperature
* Total cell energy when running in “All” mode

When “All” mode is selected, the results at each temperature point are also saved to `Results.csv`.

## Code patterns

### Temperature updates

All synthesis classes follow a common pattern for updating temperature‑dependent dictionaries:

```python
obj = Syngen(initial_T, cell)  # initialise at the first temperature
for T in temperature_range:
    obj.update_T(T)           # update thermodynamic dictionaries
    energy = obj.genome_synthesis()
```

### Dictionary structure

Thermodynamic dictionaries are keyed by metabolite code with cell‑type indexed concentration values, for example:

```python
{
    'A': {
        'name': 'dAMP',
        'MW': molecular_weight,
        'dGf[nucleotide]': gibbs_energy,
        'EC[]': ecoli_concentration,
        'SC[]': yeast_concentration,
        'M[]': mammalian_concentration,
    },
    ...
}
```

### Biological standard conversion

Conversion from chemical to biological standard is handled by subtracting the free energy of protonation at pH 7:

```python
stdbio_gibbs = std_gibbs - (N * R * T * math.log(10**-pH))
```

where `N` is the number of hydrogen atoms, `R` is the gas constant and `T` is the temperature in kelvin.

## Troubleshooting

### Reaktoro version compatibility

**Critical:** This project requires **Reaktoro v1**, not v2.  The two versions are incompatible:

* **Reaktoro v1** (1.0.7 – 1.2.3): Supports loading custom XML database files and includes organic species (amino acids, nucleic acid bases, sugars) in the SUPCRT database.
* **Reaktoro v2** (2.0+): Does not support custom XML databases and only includes inorganic species in standard SUPCRT databases.

The project depends on a custom database file (`Syncell/supcrt07-organics_AABB.xml`) that contains:
* AABB (custom amino acid backbone species)
* Organic molecules: Adenine, Cytosine, Guanine, Thymine, Uracil
* Amino acids: Alanine, Glycine, and all 20 standard amino acids
* Metabolites: Glucose, ATP, ADP, and other biological molecules

**Why not upgrade to Reaktoro v2?** Reaktoro v2 removed support for custom XML databases ([GitHub issue #222](https://github.com/reaktoro/reaktoro/issues/222)).  The standard supcrt07 database in v2 only contains ~1000 inorganic species and lacks the organic species essential for biological energy calculations.

### Common installation issues

#### "PackagesNotFoundError: reaktoro=1" on Apple Silicon

**Problem:** Reaktoro v1 was compiled for Intel Macs (x86_64) and is not available for ARM64 architecture.

**Solution:** Use the Apple Silicon installation instructions above, which run Reaktoro v1 via Rosetta 2 emulation.  This works seamlessly and has no performance issues for this application.

**Verification:**
```bash
conda activate rkt
python -c "from reaktoro import *; print('Reaktoro v1 loaded successfully')"
```

#### "RuntimeError: Could not load embedded database file"

**Problem:** You're using Reaktoro v2, which doesn't support custom XML database files.

**Solution:** Remove the v2 environment and reinstall v1:
```bash
conda deactivate
conda env remove -n rkt -y
# Then follow the installation instructions above for your platform
```

#### "Database loaded but species not found"

**Problem:** The custom database file might be corrupted or the file path is incorrect.

**Solution:** 
1. Verify the database file exists: `ls -lh Syncell/supcrt07-organics_AABB.xml`
2. Check file size (should be ~1.5 MB)
3. Run the test script: `python test_reaktoro_v1.py`
4. If issues persist, restore from the repository

#### Import errors with BioPython or other dependencies

**Problem:** Dependencies not installed or wrong Python version.

**Solution:**
```bash
conda activate rkt
pip install --upgrade biopython matplotlib progressbar2 tqdm pandas
```

### Testing your installation

A test script (`test_reaktoro_v1.py`) is provided to verify your setup.  Run it before attempting to use the main programs:

```bash
conda activate rkt
python test_reaktoro_v1.py
```

**Expected output:**
```
✓ Reaktoro v1 imported successfully
✓ Custom database loaded successfully
✓ Got Gibbs energy for AABB: -344.68 kJ/mol
✓ Got Gibbs energy for Adenine(aq): 312.84 kJ/mol
✓ Got Gibbs energy for Glycine(aq): -380.53 kJ/mol
✓ BioMolecule created for Glycine
✓ Reaktoro v1 is working correctly with your custom database!
```

If any test fails, refer to the error messages and the troubleshooting section above.

### Getting help

If you encounter issues not covered here:
1. Check that you're using Reaktoro v1 (not v2)
2. Verify your conda environment: `conda list | grep reaktoro`
3. Ensure all dependencies are installed: `pip list`
4. Review the [Reaktoro v1 documentation](https://reaktoro.org/v1/)
5. Check existing GitHub issues in the [Reaktoro repository](https://github.com/reaktoro/reaktoro/issues)
