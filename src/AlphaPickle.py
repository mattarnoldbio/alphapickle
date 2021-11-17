### AlphaPickle ###
### Version 1.3.0 ###
### Author: Matt Arnold ###
# AlphaPickle extracts results metadata from pickle (.pkl) files created by DeepMind's AlphaFold (Jumper et al., 2021, doi: 10.1038/s41586-021-03819-2)
# For detailed usage and installation instructions, please consult README.alphapickle
# New in this version: pLDDT plotting; improved plotting aesthetic

# Copyright (C) 2021  Matt Arnold

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Import dependencies
import pickle as pkl
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt ,colors as cols ,cm as cm
import json
from sys import exit

# Define class for AlphaFold metadata file and class methods
class AlphaFoldPickle:
    def __init__(self, PathToFile, FastaSequence=None, ranking=None):
        # Define attributes
        self.PathToFile = PathToFile
        self.FastaSequence = FastaSequence
        self.saving_filename = self.PathToFile.split("/")[-1].split(".")[0]
        self.saving_pathname = self.PathToFile.split(self.saving_filename)[0]
        if ranking:
            self.saving_filename = "ranked_{}".format(ranking)
        self.data = []
        self.PAE = None

        # Extract pickled data
        with (open("{}".format(self.PathToFile), "rb")) as openfile:
            while True:
                try:
                    self.data.append(pkl.load(openfile))
                except EOFError:
                    break

        # Try statement accounts for data run using non-pTM models, with no PAE output
        try:
            self.PAE = self.data[0]['predicted_aligned_error']
        except:
            print("PAE model data not present. To access this performance metric, run AlphaFold"
            "using pTM-enabled models.")

        # Define pLDDT
        self.pLDDT = self.data[0]['plddt']

    # Generate a plot from PAE measurements
    def plot_PAE(self, size_in_inches=12, axis_label_increment=100):
        ticks = np.arange(0, self.PAE[1].size, axis_label_increment)
        plt.figure(figsize=(size_in_inches,size_in_inches))
        PAE = plt.imshow(self.PAE)
        plt.xticks(ticks, fontname="Helvetica")
        plt.yticks(ticks, fontname="Helvetica")
        plt.xlabel("Residue index", size=14, fontweight="bold", fontname="Helvetica")
        plt.ylabel("Residue index", size=14, fontweight="bold", fontname="Helvetica")
        scale = plt.colorbar(PAE, shrink=0.5)
        scale.set_label(label="Predicted error (Å)",size=12, fontweight="bold", fontname="Helvetica")
        
        # Save plot
        plt.savefig('{}/{}_PAE.png'.format(self.saving_pathname, self.saving_filename),dpi=300)

        # Generate dataframe from PAE data and save to csv
        pd_PAE = pd.DataFrame(self.PAE)
        pd_PAE.to_csv('{}/{}_PAE.csv'.format(self.saving_pathname, self.saving_filename))

        return 

    # Generate a ChimeraX attribute file from pLDDT measurements   
    def write_pLDDT_file(self):
        seqMismatch = False
        pd_lDDT = pd.DataFrame(self.pLDDT)
                # Name dataframe column
        pd_lDDT.columns= ["pLDDT"]

        # If the fasta file was provided:
        if self.FastaSequence != None:
            
            # Open the fasta file in read mode
            with (open("{}".format(self.FastaSequence), "r")) as openfile:
                fasta = openfile.read()

            # Delete header line and remove line-breaks
            sequence = fasta.split("\n",1)[1].replace("\n","")

            # Check that the lengths of the two sequences match
            if len(sequence) != len(pd_lDDT):

                # If not, ignore the fasta file
                print("Length of sequence in fasta file provided ({}) does not match length of sequence used in AlphaFold prediction ({}). Ignoring fasta file.".format(len(sequence),len(pd_lDDT)))
                seqMismatch = True
            # If they do,
            else:
                # Convert the fasta sequence into a residue list
                list_sequence = []
                for item in sequence: 
                    list_sequence.append(item)

                # Convert the list into a pandas series
                pd_sequence = pd.Series(list_sequence)

                # Insert the series into the dataframe at column 1 to act as labels for the data
                pd_lDDT.insert(0,"Residue",pd_sequence)

        # Otherwise, remind user to check that they have used corret input files
        else:
            print("Number of residues for which pLDDT is provided: ", len(pd_lDDT), 
            "If this does not match the length of your sequence, please double check the input file.")



        # Tell python not to elide middle rows of dataframe when printing to std.out
        pd.set_option("display.max_rows", None, "display.max_columns", None)

        # Save dataframe to ./outputfiles with same name as original pickle and .csv extension
        pd_lDDT.to_csv('{}/{}_pLDDT.csv'.format(self.saving_pathname, self.saving_filename))
        # Delete residue ID
        if self.FastaSequence != None and seqMismatch == False:
            lDDT_table = pd_lDDT.drop('Residue', axis=1)
        else:
            lDDT_table = pd_lDDT

        # Initialise list to store Chimera-style residue identifiers (":x", where x = residue number)
        residue_list = []

        # Populate this list
        for residue in range(0,len(lDDT_table)):
            residue_list.append(":{}".format(residue+1))

        # Save to pandas format
        chimerax_numbering = pd.Series(residue_list)

        # Insert in the first column of the dataframe, to satisfy ChimeraX formatting
        lDDT_table.insert(0, 'Numbering', chimerax_numbering)

        # Tidy indices so the first label is 1 not 0
        pd_lDDT.index += 1

        # Create a file to save the Chimera attribute output into

        with (open('{}/{}_lDDT.txt'.format(self.saving_pathname, self.saving_filename),'w+')) as openfile:

            # Write file header in correct format
            openfile.write('attribute: pLDDTvalue\nmatch mode: 1-to-1\nrecipient: residues\n')

            # Iterate over rows of dataframe, writing residue ID and lDDT value to file with correct formatting
            for i,row in lDDT_table.iterrows():
                openfile.write("\t{}\t{}\n".format(row['Numbering'],row['lDDT']))

        return pd_lDDT

    # Generate a plot of pLDDT value
    def plot_pLDDT(self, size_in_inches=12, axis_label_increment=100):
        x = list(range(0,len(self.pLDDT),1))
        y = list(self.pLDDT)

        # Use standard AlphaFold colors 
        cmap = cols.LinearSegmentedColormap.from_list("", ["red","orange","yellow","cornflowerblue","blue"])

        plt.figure(figsize=(size_in_inches,(size_in_inches/2)))
        ticks = np.arange(0, len(self.pLDDT), axis_label_increment)
        plt.xticks(ticks, fontname="Helvetica")
        plt.yticks(fontname="Helvetica")
        plt.xlabel("Residue index", size=14, fontweight="bold", fontname="Helvetica")
        plt.ylabel("Predicted LDDT",size=14, fontweight="bold", fontname="Helvetica")
        plt.scatter(x,y, c=y, cmap=cmap, s=5) 
        plt.clim(0,100)
        scale = plt.colorbar(shrink=0.5)
        scale.set_label(label="Predicted LDDT",size=12, fontweight="bold", fontname="Helvetica")
        # Save to directory with pickle file in
        plt.savefig('{}/{}_pLDDT.png'.format(self.saving_pathname, self.saving_filename), dpi=300)



class AlphaFoldJson:
    def __init__(self, PathToDirectory):
        self.PathToDirectory = PathToDirectory
        self.RankingDebug = []
        try:
            with open("{}/ranking_debug.json".format(self.PathToDirectory)) as jsonfile:
                self.RankingDebugRaw = json.load(jsonfile)
            for index in enumerate(self.RankingDebugRaw['order']):
                self.RankingDebug.append(index)
        except:
            exit("To use batch processing, please ensure that the ranking_debug.json file and the result_model_n.pkl files are present in the directory issued in the command. Exiting AlphaPickle now...")

        