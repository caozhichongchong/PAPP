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
                    action='store', default=0, type=int)
parser.add_argument("-input",
                    help="The input dir or single ko list file", type=str, default='kolist',metavar='dir or file.ko')
parser.add_argument("-name",
                    help="The namefile to change the ORFs name (default: False)", type=str, default='False',metavar='filename or False')
parser.add_argument('--resultdir',
                    default="annotation", action='store', type=str, metavar='annotation',
                    help="Set the result directory for kegg identification (default: annotation)")
################################################## Definition ########################################################
args = parser.parse_args()
input_path = os.path.abspath(args.input)
if args.input_type == 1:
    in_dir, input_file = os.path.split(input_path)
    in_dir = os.path.abspath(in_dir)
    list_ko = [os.path.join(in_dir, input_file)]
else:
    in_dir = os.path.abspath(input_path)
    list_ko = glob.glob(os.path.join(in_dir, '*'))
if args.name!='False':
    name_path = os.path.abspath(args.name)
    name_dir,name_file=os.path.split(name_path)
Newname=dict()
for line in open(os.path.join(name_dir, name_file), 'rb'):
    Newname.setdefault(str(line).split('\t')[0],str(line).split('\t')[1].split('\n')[0])
try:
    os.mkdir(args.resultdir)
except OSError:
    pass
################################################### Function #########################################################
__metaclass__=type
class kegg:
    'format kegg_ko.list into table'
    def intko(self,file): ###create new method ###self is the instance itself
        f1=open('ko_formated','w')
        f1.write('KO'+'\t'+'NAME'+'\t'+'DEFINITION'+'\t'+'PATHWAY'+'\t'+'BRITE_KO1'+'\t'+'BRITE_KO2'+'\t'+'BRITE_KO3'+'\n')
        self=dict()
        Output = 0
        for line in open(file,'rb'):
            temp2 = filter(None, str(line).split(' '))
            if str(line).split(' ')[0]=='ENTRY':
                temp=temp2[1]
                self.setdefault(temp,[{'NAME':None},{'DEFINITION':None},{'PATHWAY':None},{'BRITE':None}])
            elif str(line).split(' ')[0]=='NAME':
                self[temp][0]['NAME']=' '.join(temp2[1:]).split('\n')[0]
            elif str(line).split(' ')[0]=='DEFINITION':
                self[temp][1]['DEFINITION'] = ' '.join(temp2[1:]).split('\n')[0]
            elif str(line).split(' ')[0]=='PATHWAY':
                Lable='PATHWAY'
                self[temp][2]['PATHWAY'] = [' '.join(temp2[1:]).split('\n')[0]]
            elif str(line).split(' ')[0]=='' and Lable=='PATHWAY':
                self[temp][2]['PATHWAY'].append(' '.join(temp2).split('\n')[0])
            elif str(line).split(' ')[0]=='BRITE':
                Lable='BRITE'
                IGlable=' '.join(temp2[1:]).split('\n')[0]
            elif str(line).split(' ')[0]=='' and Lable=='BRITE':
                if list(line)[11]==' ' and list(line)[12]!=' ':
                    IGlable=' '.join(temp2).split('\n')[0]
                if list(line)[12]==' ' and list(line)[13]!=' ':
                    Sublable=' '.join(temp2).split('\n')[0]
                    if 'KEGG Orthology' in IGlable: #or Transporters, Enzymes
                        self[temp][3]['BRITE']={Sublable:[]}
                elif 'KEGG Orthology' in IGlable:
                    if list(line)[13] == ' ' and list(line)[14] != ' ':
                            Subsublable = ' '.join(temp2).split('\n')[0]
                            self[temp][3]['BRITE'][Sublable].append({Subsublable: None})
                    if list(line)[14] == ' ' and list(line)[15] != ' ':
                            Subsubsublable = ' '.join(temp2).split('\n')[0]
                            for i in range(len(self[temp][3]['BRITE'][Sublable])):
                                for key in self[temp][3]['BRITE'][Sublable][i-1]:
                                    if key==Subsublable:
                                        self[temp][3]['BRITE'][Sublable][i - 1][Subsublable]=Subsubsublable
                                        f1.write(str(temp)+'\t'+str(self[temp][0]['NAME'])+'\t'+\
                                                 str(self[temp][1]['DEFINITION'])+'\t'+str(self[temp][2]['PATHWAY'])\
                                                 +'\t'+str(Sublable)+'\t'+str(Subsublable)+'\t'+str(Subsubsublable)+'\n')
                                        Output=1
            elif str(line).split(' ')[0]=='DBLINKS':
                Lable='DBLINKS'
            elif str(line).split(' ')[0]=='GENES':
                Lable='GENES'
            elif '///'in str(line):
                if Output==0:
                    f1.write(str(temp) + '\t' + str(self[temp][0]['NAME']) + '\t' + \
                         str(self[temp][1]['DEFINITION']) + '\t' + str(self[temp][2]['PATHWAY']) \
                         + '\t' + 'None' + '\t' + 'None' + '\t' + 'None' + '\n')
                Output=0
        f1.close()
    def intformatedko(self,file):
        self=dict()
        for line in open(file, 'rb'):
            temp=str(line).split('\t')[0]
            self.setdefault(temp, [{'NAME': None}, {'DEFINITION': None}, {'PATHWAY': None}, {'BRITE': None}])
            self[temp][0]['NAME'] = str(line).split('\t')[1]
            self[temp][1]['DEFINITION']= str(line).split('\t')[2]
            temp1=str(line).split('\t')[3].replace('[','').replace(']','').replace('\'','')
            self[temp][2]['PATHWAY']= temp1.split(',')
            if not self[temp][3]['BRITE']:
                self[temp][3]['BRITE'] = {str(line).split('\t')[4]: {str(line).split('\t')[5]:[str(line).split('\t')[6].split('\n')[0]]}}
            elif str(line).split('\t')[4] in self[temp][3]['BRITE']:
                if str(line).split('\t')[5] in self[temp][3]['BRITE'][str(line).split('\t')[4]]:
                    self[temp][3]['BRITE'][str(line).split('\t')[4]][str(line).split('\t')[5]].append(str(line).split('\t')[6].split('\n')[0])
                else:
                    self[temp][3]['BRITE'][str(line).split('\t')[4]].setdefault(str(line).split('\t')[5],[str(line).split('\t')[6].split('\n')[0]])
            else:
                self[temp][3]['BRITE'].setdefault(str(line).split('\t')[4],{str(line).split('\t')[5]:[str(line).split('\t')[6].split('\n')[0]]})
        return self


__metaclass__ = type
class kegg_entry:
    'get kegg'
    def getko(self,kegg_list,ko):
        self.ko=ko
        self.name=kegg_list[ko][0]['NAME']
        self.definition=kegg_list[ko][1]['DEFINITION']
        self.pathway=kegg_list[ko][2]['PATHWAY']
        self.brite = kegg_list[ko][3]['BRITE']


################################################### Programme #########################################################
###format kegg file
filename='ko'
Kegg=kegg()
Kegg.intko(filename)

###read ko_formated
filename='ko_formated'
Kegg=kegg()
Kegg=Kegg.intformatedko(filename)
for file_name in list_ko:
    try:
        in_dir, input_file = os.path.split(file_name)
        f2 = open(os.path.join(args.resultdir, os.path.splitext(str(input_file))[0] + '.kegg'), 'w')
        f2.write('Genome\tORF\tKO\tNAME\tDEFINITION\tPATHWAY\tBRITE_KO1\tBRITE_KO2\tBRITE_KO3\n')
        KO=kegg_entry()
        KO = kegg_entry()
        for line in open(os.path.join(in_dir, file_name), 'rb'):
            if len(str(line).split('\t'))>1:
                if args.name != 'False':
                    for oldname in Newname:
                        if oldname in str(line).split('\t')[0]:
                            newname = str(line).split('\t')[0].replace(oldname, Newname.get(oldname))
                else:
                    newname=str(line).split('\t')[0]
                try:
                    KO.getko(Kegg, str(line).split('\t')[1].split('\n')[0])
                    for Sublable in KO.brite:
                        for Subsublable in KO.brite[Sublable]:
                            f2.write(str(newname).split('_')[0].split('-')[0] + '\t' + str(newname) + '\t' + KO.ko + \
                                     '\t' + KO.name + '\t' + KO.definition + '\t' + str(KO.pathway)+'\t'+str(Sublable) + '\t' + \
                                     str(Subsublable) + '\t' + str(KO.brite[Sublable][Subsublable]) + '\n')
                except KeyError:
                    print str(line).split('\t')[1].split('\n')[0]+'missing!'
    except IOError:
        print 'Files were missing?'
    finally:
        print 'Cleaning Up...'
        del file_name
