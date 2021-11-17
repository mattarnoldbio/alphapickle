# AlphaPickle
Version : 1.4.0
Author : Matt Arnold (m.arnold.1@research.gla.ac.uk)
Date : 16-Nov-21

ALPHAPICKLE is a programme for extracting the outputs of DeepMind's ALPHAFOLD protein prediction algorithm (Jumper et al., 2021, 10.1038/s41586-021-03819-2) from the pickle (.pkl) file output.
- Currently, the outputs of the programme are: a plot of the predicted aligned error (PAE), a ChimeraX attribute file for colouring the outputted models by pLDDT value and csv files containing the predicted aligned error and per-residue confidence in the predicted model (expressed as local distance difference test (lDDT) values, as well as plots of these two metrics, saved as .png files; see Mariani et al., 2013 10.1093/bioinformatics/btt473). Consult the original publication for details on the importance of these data. Other metadata contained in the .pkl are described at  https://github.com/deepmind/alphafold/blob/main/README.md. If you would like to see any of these added to future version, please get in contact.
- For details on how to setup and use ALPHAPICKLE, see below.

# Use

- The program has two basic options: input a single metadata file (by issuing the the absolute path to the file with the flag -pf), or process metadata for all models output by a run of AlphaFold (by issuing the absolute path to the results directory, which must contain both the ranking_debug.json file and all of the result_model*.pkl files). For each metadata file processed, the outputs (plots and tables of pLDDT and PAE, as well as a ChimeraX attribute file for pLDDT) are saved to the same directory as the input files.
- In addition to the basic options, advanced options for tinkering with the dimensions/scaling of plots are provided. Details are available on issuing run_AlphaPickle.py -h. If the plot is not sufficient after adjusting these parameters, the user may find it easier to recreate the plot from the outputted csv files. 

# Install

- Python dependencies are provided in the requirements.txt file. 

# Citation

- If you use AlphaPickle in your work (during analysis, or for plots that end up in publications), please cite AlphaPickle as follows: Arnold, M. J. (2021) AlphaPickle doi.org/10.5281/zenodo.5708709

[![DOI](https://zenodo.org/badge/429171188.svg)](https://zenodo.org/badge/latestdoi/429171188)
