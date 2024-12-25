# Environmental transitions and HGT rates

Code to reproduce the results of the paper, "Horizontal gene transfer rates decrease during transitions to new ecosystems".

Jupyter notebooks in `notebooks` directory were run in the order of their numbering. Scripts and other helper functions that they depend on or are mentioned in the notebooks can be found in the `code` directory.

Please note that using 'Run all' or equivalent in the Jupyter notebooks will generally not be useful. Some of the intervening steps in the notebooks are markdown cells instructing how to run programs via shell, separately. These programs are time consuming ones, often utilising multiprocessing in an HPC. Please run the notebooks, cell by cell, keeping this in mind. Generally speaking, if you have run the previuos cell, you can run the current cell if it's a Python cell. If it's a markdown cell instructing you to do something, do it before proceeding.

This study downloads and processes an number of large files, which were stored in the `data` directory, which is empty in this repo. If you follow/run the notebooks you will progressively fill the `data` directory to reproduce all the results as well a the figures in the paper. Alternatively, you can download this repo containing the results and figures in the `data` directory, from the Zenodo link in the paper.

## Required python packages

A Mamba/Conda environment called `hgt_analyses` was used for all the analyses. This environment with all required packages can easily be created again using the `mamba_packages.yml` file, using the following command:
```
mamba env create -f mamba_packages.yml
```
I use Mamba since it's faster than using regular Conda but if you use Conda just replace `mamba` with `conda` in the above command.

The only other file you need to install is `xlsx2csv` using `pip`:
```
pip install xlsx2csv
```

Additionally, I used `ripgrep` instead of `grep` for faster extraction of subset data from EggNOG. In case you rely on `grep` make edits accordingly in the notebooks (replace commands of `rg`)
