# Panpath
#Copyright: LG209, Environmental Biotechnology Lab, The University of Hong Kong
#Author: An Ni ZHANG, email: caozhichongchong@gmail.com
#Panpath for new Accumulibacter (Acc) genomes and new Acc transcriptional datasets visualization 
#in Acc reference pan-pathway, specifically for anaerobic (AN) and aerobic (AE) phases;

#Supplementary files
#Edge_number. Empty edge table
#Node_number. Empty node table
#File S3. Reference pan-pathway of Acc in AN phase for LCRPKM values
#File S4. Reference pan-pathway of Acc in AE phase for LCRPKM values
#File S5. Reference pan-pathway of Acc in AN phase for CRPKM ranks
#File S6. Reference pan-pathway of Acc in AE phase for CRPKM ranks

#1. Files S3 and S4 for CRPKM values of metatranscriptomid datasets.
#Files S5 and S6 for CRPKM ranks of metatranscriptomid datasets
#2. Before any changes, make sure the File S1 (empty edge table) and File S2 (empty node table) are imported
#3. running PAO_Panpathway.py, you may firstly try with the Example.txt.
python PAO_Panpathway.py -i Example.txt -ko 1 -col 2,3

usage: PAO_Panpathway.py [-h]
                         [-i your file with KO numbers and their color lable integer]
                         [-ko 1] [-col 2,3,4]
                         [--f tab: default or $'	' , blackspace: ' ']
                         [--p 0 or 1 or 2] [--v r or v]

optional arguments:
  -h, --help            show this help message and exit
  -i your file with KO numbers and their color lable (integer)
                        A single file name
  -ko 1                 The colomn number of KO
  -col 2,3,4            The colomns number of expression levels (separated by
                        ',')
  --f tab: default or $'	' , blackspace: ' '
                        The seperation character of input file (default: tab)
  --p 0 or 1 or 2       The method for merging labels for path with multiple
                        subpaths (0 for common, 1 for max, 2 for min)
  --v r or v            CRPKM ranks or CRPKM value (default)(r for CRPKM
                        ranks, v for CRPKM value)
#The CRPKM value would be transformed by log10
#4. After running PAO_Panpathway.py, open Cytoscape (File S3-S6)
#Please import the edge (to the Edge Table Colomn) and node table (to the Node Table Colomn) of your datasets
#The Key Colomn for Network must be set to 'No'
#The attribute of the colomn 'Color' should be set to 'string', not 'int' or 'float'
