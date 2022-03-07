### AlphaPickle ###
### Version 1.2.1 ###
### Author: Matt Arnold ###
# AlphaPickle extracts results metadata from pickle (.pkl) files created by DeepMind's AlphaFold (Jumper et al., 2021, doi: 10.1038/s41586-021-03819-2)
# For detailed usage and installation instructions, please consult README.alphapickle
# New in this version: Bug fix @ line43: Conditional was non-functional. Replaced with working version.
#                       Bug fix @ line39: incomplete {}.format construction led to no function. remedied.
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

import src.AlphaPickle as ap
import argparse
import glob
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="AlphaPickle\n"
    "Version 1.4.0\n"
    "Input the location of an AlphaFold output directory and generate \n"
    "a PAE plot (if pTM models were used), a pLDDT plot and a ChimeraX attribute file \n"
    "containing pLDDT data for each model. Both of these metrics are also exported to \n"
    "csv files. All outputs save to the  directory containing the input files.\n"
    "Copyright (C) 2021  Matt Arnold")
    parser.add_argument("-od","--output_directory", help='Path to AlphaFold output directory')
    parser.add_argument("-ps","--plot_size", help='Optional (Default = 12). Change size (in inches) of plots. This may be useful for very short or long input sequences', default=12)
    parser.add_argument("-pi","--plot_increment", help='Optional (Default = 100). Change the increment of plot axis labels using residue numbering. This may be useful for very short or long input sequences', default=100)
    args = parser.parse_args()
    
    print("Launching AlphaPickle!")
    print("\n")

    print("""
        _       _             ____ ___ ____ _  ___     _____ 
   __ _| |_ __ | |__   __ _  |  _ \_ _/ ___| |/ / |   | ____|
  / _` | | '_ \| '_ \ / _` | | |_) | | |   | ' /| |   |  _|  
 | (_| | | |_) | | | | (_| | |  __/| | |___| . \| |___| |___ 
  \__,_|_| .__/|_| |_|\__,_| |_|  |___\____|_|\_\_____|_____|
         |_|                                                 
                                 
    """)

    rankings = ap.AlphaFoldJson(args.output_directory).RankingDebug
    for model in rankings:
        print("Processing ranked model {} (result_{}).".format(str(model[0]),model[1]))
        results = ap.AlphaFoldPickle(args.output_directory + "/result_" + model[1] + ".pkl",args.fasta_file,ranking=str(model[0]))
        results.write_pLDDT_file()
        results.plot_pLDDT(size_in_inches=args.plot_size, axis_label_increment=args.plot_increment)
        if type(results.PAE) == np.ndarray:
            results.plot_PAE(size_in_inches=args.plot_size, axis_label_increment=args.plot_increment)
        print("\n")
        
    print("Processing complete!")
    print("\n")
    print("Data saved to output directory")