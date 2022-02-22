# AlphaPickle [![DOI](https://zenodo.org/badge/429171188.svg)](https://zenodo.org/badge/latestdoi/429171188)
Version : 1.5.0

Author : Matt Arnold (m.arnold.1@research.gla.ac.uk)
Date : 16-Nov-21

ALPHAPICKLE is a programme for extracting the outputs of DeepMind's ALPHAFOLD protein prediction algorithm (Jumper et al., 2021, 10.1038/s41586-021-03819-2) from the pickle (.pkl) file output.
- Currently, the outputs of the programme are: a plot of the predicted aligned error (PAE), a ChimeraX attribute file for colouring the outputted models by pLDDT value and csv files containing the predicted aligned error and per-residue confidence in the predicted model (expressed as local distance difference test (lDDT) values, as well as plots of these two metrics, saved as .png files; see Mariani et al., 2013 10.1093/bioinformatics/btt473). Consult the original publication for details on the importance of these data. Other metadata contained in the .pkl are described at  https://github.com/deepmind/alphafold/blob/main/README.md. If you would like to see any of these added to future version, please get in contact.
- For details on how to setup and use ALPHAPICKLE, see below.

# Use

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
 - Bug fix at /src/AlphaPickle.py-line163: fix broken reference to dataframe.
