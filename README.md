# AlphaPickle [![DOI](https://zenodo.org/badge/429171188.svg)](https://zenodo.org/badge/latestdoi/429171188)
Version : 1.5.1

Author : Matt Arnold (m.arnold.1@research.gla.ac.uk)
Date : 16-Nov-21

ALPHAPICKLE is a programme for extracting the outputs of DeepMind's ALPHAFOLD protein prediction algorithm (Jumper et al., 2021, 10.1038/s41586-021-03819-2).

AlphaPickle is multipurpose Python script for producing plots and user-legible files from the output of AlphaFold2 (notebook) and Colabfold (notebook).

The functions provided here aim to tackle the problem of output metadata from these programs being difficult to process for users who don't have Python training. Currently, PAE outputs from AlphaFold 2.1.0 are in the form of .pkl files, which can only be read using a Python script. For data from AlphaFold DB and ColabFold, these data are in .json format (not typically easy to process for non-code writing users).

For all of the above software, pLDDT values are outputted in the B-factor field of the PDB file for each prediction. This is very useful for visualisation (e.g. using the ChimeraX command color bfactor palette alphafold), but may be difficult in terms of customisable visualisation for non-code writing users.

This collection of code will take any of the above output files and provide a .csv file (which can be opened and used for plotting in Excel, Numbers, Google Sheets) as well as a downloadable plot.

- Currently, the outputs of the programme are: csv files containing the predicted aligned error and per-residue confidence in the predicted model (expressed as local distance difference test (lDDT) values, as well as plots of these two metrics, saved as .png files; see Mariani et al., 2013 10.1093/bioinformatics/btt473). Consult the original publication for details on the importance of these data. Other metadata contained in the .pkl are described at  https://github.com/deepmind/alphafold/blob/main/README.md. If you would like to see any of these added to future version, please get in contact.
- For details on how to setup and use ALPHAPICKLE, see below.

# Use

- EASIEST: visit the Colab page at https://colab.research.google.com/github/mattarnoldbio/alphapickle/blob/main/AlphaPickle.ipynb
    - No command line interaction required
    - No system configuration or Python package installation required
    - Local installation is now probably only necessary if only .pkl files are available for metadata, as the size of these may be prohibitively large when considering upload to Colab.

- Local use:
- The program has two basic options: input a single metadata file (by issuing the the absolute path to the file with the flag -pf), or process metadata for all models output by a run of AlphaFold (by issuing the absolute path to the results directory, which must contain both the ranking_debug.json file and all of the result_model*.pkl files). For each metadata file processed, the outputs (plots and tables of pLDDT and PAE, as well as a ChimeraX attribute file for pLDDT) are saved to the same directory as the input files.
- In version 1.5.0, functionality has been added to allow processing of the outputs of DeepMind's Colab notebook (plotting pLDDT from the b-factor column of an AlphaFold-generated PDB file and PAE from a JSON file). This function will also be useful if using data from AlphaFoldDB.
- Usage examples:
    - To process all metadata files in an AlphaFold results directory (recommended; requires that directory also contains raking_debug.json file): `python3 run_AlphaPickle.py -od /absolute/path/to/output/directory`
    - To process a specific file: `python3 run_AlphaPickle.py -pf /absolute/path/to/pickle/file`
    - To produce a pLDDT plot from an AlphaFold PDB file: `python3 run_AlphaPickle.py -pdb /absolute/path/to/pdb/file`
    - To produce a PAE plot from an AlphaFold Colab or AlphaFold DB .json file: `python3 run_AlphaPickle.py -json /absolute/path/to/predicted_aligned_error.json/file`
    - n.b. If `python3` does not point to your system's binary executable file for Python 3, replace `python3` in the above commands with this. This software has only been tested with Python 3.9.5 and support for other versions of Python is not guaranteed. If you do not add the location of `run_AlphaPickle.py` to your system's `$PATH` variable, you will need to run these scripts from a terminal in the directory containing `run_AlphaPickle.py`

- In addition to the basic options, advanced options for tinkering with the dimensions/scaling of plots are provided. Details are available on issuing run_AlphaPickle.py -h. If the plot is not sufficient after adjusting these parameters, the user may find it easier to recreate the plot from the outputted csv files. 

# Install

- How you install Python dependencies on your system is obviously up to you, but the relevant modules are listed in *requirements.txt files. Currently, this program has only been tested on OSX and Linux (if you test on Windows and find problems, please raise an issue). The module dependencies are subtly different between OS's, thanks to the long list of dependencies for matplotlib. Please make sure you select the *requirements.txt file relevant to your OS (either osx64_requirements.txt for MacOS or linux_requirements.txt for Linux distros). Consider installing the dependencies in a virtual environment to avoid potential version conflicts.
 

# Citation

- If you use AlphaPickle in your work (during analysis, or for plots that end up in publications), please cite AlphaPickle as follows: Arnold, M. J. (2021) AlphaPickle doi.org/10.5281/zenodo.5708709

# Changes
 - pLDDT to csv from pdb file
