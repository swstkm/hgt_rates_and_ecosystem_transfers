# Environmental transitions and HGT rates

Code to reproduce the results of the paper, "Streamlined genomes, not horizontal gene transfer, mark bacterial transitions to unfamiliar environments" (BioRxiv preprint [here](https://www.biorxiv.org/content/10.1101/2024.12.27.630308)).

This repository is only for analyses of ecosystem transfers and associated HGT rates. For the corresponding analysis about pathogencity gains instead, see [this repo instead](https://gitlab.cs.uni-duesseldorf.de/mishra/hgt_rates_and_pathogenicity_gains)

Jupyter notebooks in `notebooks` directory were run in the order of their numbering. Scripts and other helper functions that they depend on or are mentioned in the notebooks can be found in the `code` directory.

Please note that using 'Run all' or equivalent in the Jupyter notebooks will generally not be useful. Some of the intervening steps in the notebooks are markdown cells instructing how to run programs via shell, separately. These programs are time consuming ones, often utilising multiprocessing in an HPC. Please run the notebooks, cell by cell, keeping this in mind.

This study downloads and processes an number of large files, which were stored in the `data` directory, which is empty in this repo. If you follow/run the notebooks you will progressively fill the `data` directory to reproduce all the results as well a the figures in the paper. Alternatively, you can download this repo containing the results and figures in the `data` directory, from the Zenodo link in the paper.

## Required python packages

A Mamba/Conda environment called `hgt_analyses` was used for all the analyses. This environment with all required packages can easily be created again using the `mamba_packages.yml` file, using the following command:
```
mamba env create -f mamba_packages.yml
```

The only other file you need to install is `xlsx2csv` using `pip`:
```
pip install xlsx2csv
```

Additionally, I used `ripgrep` instead of `grep` for faster extraction of subset data from EggNOG. In case you rely on `grep` make edits accordingly in the notebooks (replace commands of `rg`)
