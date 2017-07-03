import os
from Bio import SeqIO
import glob
import argparse

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-input_type",
                    help="The input type, eg: 0 for a dir; 1 for a file",
                    metavar="0 or 1",
                    choices=[0, 1],
                    action='store', default=1, type=int)
parser.add_argument("-input",
                    help="The input dir or a single blast file", type=str, default='blast.txt',metavar='dir or blast.txt')
parser.add_argument('--resultdir',
                    default="kolist", action='store', type=str, metavar='kolist',
                    help="Set the result directory for kegg identification (default: ko)")
parser.add_argument('--blast_iden',
                    default=80, action='store', type=float, metavar='80',
                    help='Identity cutoff for kegg annotation (default is 80)')
parser.add_argument('--blast_leng',
                    default=0.8, action='store', type=float, metavar='0.8',
                    help='Hitlength cutoff for kegg annotation (default is 0.8)')
################################################## Definition ########################################################
args = parser.parse_args()
input_path = os.path.abspath(args.input)
if args.input_type == 1:
    in_dir, input_file = os.path.split(input_path)
    in_dir = os.path.abspath(in_dir)
    list_blast = [os.path.join(in_dir, input_file)]
else:
    in_dir = os.path.abspath(input_path)
    list_blast = glob.glob(os.path.join(in_dir, '*'))
try:
    os.mkdir(args.resultdir)
except OSError:
    pass
Cutoff_identity = float(args.blast_iden)
Cutoff_hitlength = float(args.blast_leng)
Length=dict()
for line in open('prokaryotes_pep_length.txt','rb'):
    Length.setdefault(str(line).split('\t')[0],int(str(line).split('\t')[2].split('\n')[0]))
KO=dict()
for line in open('ko_genes.list','rb'):
    KO.setdefault(str(line).split('\t')[1].split('\n')[0],str(line).split('\t')[0].split(':')[1])
try:
    os.mkdir(args.resultdir)
except OSError:
    pass
################################################### Programme #########################################################
for file_name in list_blast:
    try:
        in_dir, input_file = os.path.split(file_name)
        f1 = open(os.path.join(args.resultdir, os.path.splitext(str(input_file))[0] + '.ko'), 'w')
        for line in open(os.path.join(in_dir, input_file), 'rb'):
            if float(str(line).split('\t')[2])>=Cutoff_identity:
                if float(str(line).split('\t')[3])/float(Length[str(line).split('\t')[1]])>=Cutoff_hitlength:
                    try:
                        f1.write(str(line).split('\t')[0]+'\t'+KO[str(line).split('\t')[1]]+'\n')
                    except KeyError:
                        print str(line).split('\t')[1]
        f1.close()
    except IOError:
        print 'Files were missing?'
    finally:
        print 'Cleaning Up...'
        del file_name


