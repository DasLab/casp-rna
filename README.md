![Alt text](img/logo.png "Title")

#

**casp-rna is the official CASP15 pipeline to assess the accuracy of submitted models for RNA structure prediction.**



# About

This repository contains the code for the casp-rna pipeline, a tool developed to assess the accuracy of submitted models for RNA structure prediction. The pipeline calculates ZRNA, a weighted Z-score average of several different assessment metrics, to evaluate the models. To capture the topology, local environment, and geometries of RNA, ZRNA incorporates additional metrics beyond RMSD, including TM-score, GDT-TS, INF scores, lDDT, and clashscore. The pipeline encompasses a workflow for data wrangling, job parallelization, and ranking visualizations. A two-pass procedure is employed for Z-scores, and models with initial Z-scores falling under a tolerance threshold of -2 are discarded in the first pass. The pipeline also compares submitted predictor models to all available experimental models for RNA with multiple conformations.

# Organization

The CASP-RNA repository contains all scripts and code used to obtain and analyse scores for the CASP15 RNA category. The repository is organised as follows:
<!-- 
- `bins` contains the binaries used in the analysis, or a symlink to the binaries (if those packages are already installed on the system).
- `downloads` contains
- `scripts` contains the scripts used to run the analysis.
- `scores` contains results and summary of Z-rna and other metrics.
-  -->

# Installation

Clone the repository to your local machine:

```
$ git clone https://github.com/philpham8/casp-pipeline.git
$ cd casp-rna
```

## General dependencies:
- `python >= 3.4`
- `py-numpy/1.20.3_py39`
- `scipy`
- `py-matplotlib/3.4.2_py39`
- `py-pandas/1.3.1_py39`
- `seaborn`
- `biopython`
- `viz`
- `py-scipy/1.6.3_py39`

## Dependencies for metrics calculations:
- `US-align` for TM-score calculations
- `LGA` for GDT-TS calculations
- `PHENIX` for molprobity clashscore calculations
- `OpenStructure` for lDDT calcuations
- `rna-tools` for ClaRNA-based INF scores calculations

Streamline installation of the five packages required for metrics calculations is provided in the `setup.sh` script. The script will install the required packages to the `bins/` directory. For LGA and PHENIX, their binaries must be obtained through their website and placed into casp-rna's `downloads/` folder. Subsequently, the script can be run as follows:

```
$ bash setup.sh
```

Existing installations of the above packages can be used. When running `setup.sh` script, the script can create a symlink to the existing installation of the packages. Alternatively, the script can install the packages in the `bins/` directory.


### PHENIX 

Installation of PHENIX is described in the [PHENIX documentation](https://www.phenix-online.org/documentation/reference/installation.html). PHENIX can be downloaded [here](https://phenix-online.org/download/) and placed in the project's `downloads/` folder. 

### LGA

Binary for LGA can be downloaded [here](http://proteinmodel.org/AS2TS/LGA/lga.html). The binary must be placed in the project's `downloads/` folder before running the `setup.sh` script.

# Usage

The scripts can be run from the command line, or imported into a python script. The scripts can be run from the command line, or imported into a python script.

### Structure of input directory
Each target directory must contain two subdirectories: `references/` (ground truth models) and `models/` (predicted models). The `references/` and `models/` subdirectories each must contain at least one .pdb file. Different states and configurations can be stored in separate files. The target directory can be placed anywhere on the system, but the path to the target directory must be provided to the script.

### Output
Metric summary is exported to a .csv file in the `scores/` directory, and graphs are exported to the `figures/` directory. Intermediate files for debugging for further data exploration are stored in the `runs/` directory.

### `parallel.py`

This script is used to run the pipeline on a computing cluster. The script takes in a list of target directories and runs the pipeline on each target directory. The script can be run as follows:

```
bash parallel/parallel.sh {path_to_target} {metric}
```

where {path_to_target} is the relative path containing the `references/` and `models/` subdirectories and {metric} is the desired metric to be calculated ("inf", "gdt", "tm_score", "lddt", "clashscore", "rmsd", or "all").

Execution of `parallel.sh` will launch a large number of jobs. Some minor modifications to match to the specification of your institution's computing cluster job scheduler may be required.

### `make_figures.ipynb`

This notebook contains the code to generate the figures in the paper.



# Reference

If you use this code, please cite the following paper:

    @article{casp15_rna,
        author = {TBD},
        title = {Assessment of three-dimensional RNA structure prediction in CASP15},
        journal = {TBD},
        year = {TBD},
        volume = {TBD},
        number = {TBD},
        pages = {TBD},
        doi = {TBD},
        url = {TBD}
    }