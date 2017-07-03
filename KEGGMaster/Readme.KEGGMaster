# KEGGMaster
#Copyright: LG209, Environmental Biotechnology Lab, The University of Hong Kong
#Author: An Ni ZHANG, email: caozhichongchong@gmail.com
#KEGGMaster consists of three python scripts to assign KO numbers to sequences and visualize metabolism and non-mtabolism pathways.
#The combination of three python scripts will convert blast result into highlighted kgml files
#The combination of Kegg_Ano.py and Kegg_Parse.py will convert KOlist into highlighted kgml files

#1.Blast result filtering (Kegg_blast.py) generate KOlist based on the blast result in txt format, a length file of KO sequences
#-input_type 0 for a dir, 1 for a file; -input yourdir or file; --blasst_iden default 80; --blast_leng default 0.8
python Kegg_blast.py -input_type 0 -input inputdir

#2.KO annotation (Kegg_Ano.py) formated the ko file (from Kegg Database) and assign pathway information to each KO entry
#ko_map and pathway.list are requred from Kegg Database
#-input_type 0 for a dir, 1 for a file; -input yourdir or file; -name a file with old and new names to change ORF names, optional;
python Kegg_Ano.py -input_type 0 -input yourdir -name Genome_change_name.txt

#3.KGML pathway generation (Kegg_Parse.py) to convert KOlist into requirint ko_gene.list and ko_list frm Kegg Database
#3-0. python Kegg_Parse.py -anno_type 0 
#This step generates yourfile.overall.pathway, yourfile.ko.pathway (KO in each pathway), yourfile.ko.gene.pathway (KO and genes in each pathway)
#-input_type 0 for a dir, 1 for a file; -input yourdir or file; --gene input colomn of ORFs; --ko input column of KO; --format the separation character of your input file: "\t" or " "
python Kegg_Parse.py -input_type 0 -input yourdir -anno_type 0 --gene ORFcolumn --ko KOcolumn --format " "

#3-1. python Kegg_Parse.py  -anno_type 1
#You may edit the overall.pathway to simply delete the rows of pathways you donnot with to display in advance
#kgml files containing all metabolism and non-metabolism pathways from Kegg Database are required
python Kegg_Parse.py -input_type 0 -input yourdir -anno_type 1

#4. KGML visualization
#KEGGScape in Cytoscape required
#yourfile.ko00110.meta and yourfile.ko00110.meta.filter for overall metabolism pathways
#yourfile.koXXXXX.meta for metabolism pathways, yourfile.koXXXXX.other for non-metabolism pathways
Cytoscape > network import > Kegg database
