import os
from Bio import SeqIO
import glob
import argparse

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-anno_type',
                    help="The annotation type, eg: 0 for pathway annotation of ko_list; 1 for kgml generation of path_list (default: 0)",
                    metavar="0 or 1",
                    choices=[0, 1],
                    action='store', default=0, type=int)
parser.add_argument("-input_type",
                    help="The input type, eg: 0 for a dir; 1 for a file",
                    metavar="0 or 1",
                    choices=[0, 1],
                    action='store', default=1, type=int)
parser.add_argument("-input",
                    help="The input dir or single file", type=str, default='ko_list',metavar='dir or file')
parser.add_argument("--gene",
                    help="The colomn number of gene names in input ko_list", type=int, default=1,metavar='1')
parser.add_argument("--ko",
                    help="The colomn number of ko in input ko_list", type=int, default=2,metavar='2')
parser.add_argument("--format",
                    help="The seperation character of input ko_list (default: tab)", type=str, default="\t",metavar='tab: "\\t", blackspace: " "')
parser.add_argument('--keggdir',
                    default="kegg", action='store', type=str, metavar='kegg',
                    help="Set the directory for kegg pathway list (default: kegg)")
parser.add_argument('--kgmldir',
                    default="kgml", action='store', type=str, metavar='kgml',
                    help="Set the directory for kgml of kegg pathway (default: kgml)")
parser.add_argument('--resultdir',
                    default="result", action='store', type=str, metavar='result',
                    help="Set the result directory for integron identification (default: result)")
################################################## Definition ########################################################
args = parser.parse_args()
input_path = os.path.abspath(args.input)
kegg_path = os.path.abspath(args.keggdir)
if args.input_type == 1:
    in_dir, input_file = os.path.split(input_path)
    in_dir = os.path.abspath(in_dir)
    list_file = [os.path.join(in_dir, input_file)]
else:
    in_dir = os.path.abspath(input_path)
    list_file = glob.glob(os.path.join(in_dir,'*'))
try:
    os.mkdir(args.resultdir)
except OSError:
    pass
Anno_type=int(args.anno_type)

################################################### Function ########################################################
def intkolist(file,Col_gene,Col_ko,Sept):
    KO=dict()
    try:
        for line in open(os.path.join(in_dir, str(file)), 'rb'):
            if KO.get(str(line).split(Sept)[Col_ko].split('\n')[0].split('\r')[0].replace('\"',''),'None')=='None':
                KO.setdefault(str(line).split(Sept)[Col_ko].split('\n')[0].split('\r')[0].replace('\"',''),
                              [str(line).split(Sept)[Col_gene].split('\n')[0].split('\r')[0]])
            else:
                KO[str(line).split(Sept)[Col_ko].split('\n')[0].split('\r')[0].replace('\"','')].append(str(line).split(Sept)[Col_gene].split('\n')[0].split('\r')[0])
        return KO
    except IOError:
        print 'Files were missing?'


def pathassign(KO,Path):
    Pathway=dict()
    for ko in KO:
        if Path.get(ko,'None')!='None':
            for path in Path[ko]:
                if Pathway.get(path,'None')=='None':
                    Pathway.setdefault(path,[ko])
                else:
                    Pathway[path].append(ko)
        elif Pathway.get('None','None')=='None':
            Pathway.setdefault('None',[ko])
        else:
            Pathway['None'].append(ko)
    return Pathway


def writepath(Pathway,KO):
    Dir, File = os.path.split(file_name)
    f1 = open(os.path.join(args.resultdir, os.path.splitext(str(File))[0] + '.overall.pathway'), 'wb')
    f2 = open(os.path.join(args.resultdir, os.path.splitext(str(File))[0] + '.ko.pathway'), 'wb')
    f3 = open(os.path.join(args.resultdir, os.path.splitext(str(File))[0] + '.ko.gene.pathway'), 'wb')
    for line in open(os.path.join(kegg_path, 'pathway.list'), 'rb'):
        if '#' in str(line):
            f1.write(str(line))
            f2.write(str(line))
            f3.write(str(line))
        elif str(line).split('\t')[0] in Pathway:
            temp=Pathway[str(line).split('\t')[0]]
            f1.write(str(line).replace('\n', '\t') + '(' + str(len(temp)) + ')\n')
            for ko in temp:
                f2.write(str(line).replace('\n','\t')+'('+str(len(temp))+')\t'+str(ko)+'\n')
                f3.write(str(line).replace('\n', '\t') + '(' + str(len(temp)) + ')\t' + str(ko) +'\t'+str(KO[ko])+ '\n')
    f1.write('#None\n')
    f2.write('#None\n')
    f3.write('#None\n')
    temp = Pathway['None']
    f1.write('None' + '(' + str(len(temp)) + ')\n')
    for ko in temp:
        f2.write('None' + '(' + str(len(temp)) + ')\t' + str(ko) + '\n')
        f3.write('None' + '(' + str(len(temp)) + ')\t' + str(ko) + '\t' + str(KO[ko]) + '\n')
    f1.close()
    f2.close()
    f3.close()


def pathlist(file_name):
    try:
        KO_meta=[]
        KO_other=[]
        Path_meta=[]
        Path_other=[]
        Dir, File = os.path.split(file_name)
        for line in open(os.path.join(args.resultdir, str(os.path.splitext(str(File))[0]))+'.overall.pathway', 'rb'):
            if "#" in str(line) and "##"not in str(line):
                if 'Metabolism' in str(line):
                    Lable='Meta'
                else:
                    Lable='Other'
            elif "#" not in str(line):
                if Lable=='Other':
                    Path_other.append(str(line).split('\t')[0])
                else:
                    Path_meta.append(str(line).split('\t')[0])
        for line in open(os.path.join(args.resultdir, str(os.path.splitext(str(File))[0]))+'.ko.pathway', 'rb'):
            if "#" in str(line) and "##"not in str(line):
                if 'Metabolism' in str(line):
                    Lable='Meta'
                else:
                    Lable='Other'
            elif "#" not in str(line):
                if Lable=='Other':
                    if str(line).split('\t')[0] in Path_other:
                        KO_other.append(str(line).split('\t')[-1].split('\r')[0].split('\n')[0])
                else:
                    if str(line).split('\t')[0] in Path_meta:
                        KO_meta.append(str(line).split('\t')[-1].split('\r')[0].split('\n')[0])
        for path in Path_other:
            f1=open(os.path.join(args.resultdir, str(os.path.splitext(str(File))[0])+'.ko'+str(path)+'.other.xml'), 'wb')
            Newentry=0
            for line in open(os.path.join(args.kgmldir, 'ko'+str(path)+'.xml'), 'rb'):
                if '</' in str(line):
                    Newentry=0
                    f1.write(str(line))
                elif 'name="ko:' in str(line):
                    if any(x in str(line).split('name="')[1].split('"')[0] for x in KO_other):
                        Newentry=1
                    f1.write(str(line))
                elif "bgcolor" in str(line) and Newentry==1:
                    f1.write(str(line).split('fgcolor="')[0] + 'fgcolor="#FF0000" bgcolor="#FF0000"' +
                             str(line).split('bgcolor="')[1].split('"')[1])
                else:
                    f1.write(str(line))
            f1.close()
        Cpd_meta = []
        Reac_meta=[]
        for path in Path_meta:
            f1 = open(os.path.join(args.resultdir, str(os.path.splitext(str(File))[0])+'.ko' + str(path) + '.meta.xml'), 'wb')
            Newentry = 0
            for line in open(os.path.join(args.kgmldir, 'ko'+str(path)+'.xml'), 'rb'):
                if '</' in str(line):
                    Newentry=0
                    f1.write(str(line))
                elif 'name="ko:' in str(line):
                    if any(x in str(line).split('name="')[1].split('"')[0] for x in KO_meta):
                        Newentry=1
                        try:
                            temp=str(line).split('reaction="')[1].split('"')[0]
                            for Reaction in temp.split(' '):
                                Reac_meta.append(Reaction)
                        except IndexError:
                            pass
                    f1.write(str(line))
                elif 'name="cpd:' in str(line) or 'name="gl:' in str(line):
                    Cpd_meta.append(str(line).split(':')[1].split('"')[0])
                    f1.write(str(line))
                elif "bgcolor" in str(line) and Newentry==1:
                    f1.write(str(line).split('fgcolor="')[0] + 'fgcolor="#FF0000" bgcolor="#FF0000"'+
                             str(line).split('bgcolor="')[1].split('"')[1])
                else:
                    f1.write(str(line))
            f1.close()
        f1 = open(os.path.join(args.resultdir, str(os.path.splitext(str(File))[0])+'.ko01100.meta.xml'), 'wb')
        Newentry = 1
        Newtag=1
        Newcpd=1
        Reac_meta_filter=[]
        Cpd_meta_filter=[]
        for line in open(os.path.join(args.kgmldir, 'ko01100.xml'), 'rb'):
            if ('name="cpd:' in str(line) or 'name="gl:' in str(line)) and 'entry' in str(line):
                Newentry=1
                f1.write(str(line))
                if str(line).split(':')[1].split('"')[0] in Cpd_meta:
                    Newtag = 1
                    Lable='cpd'
                else:
                    Newtag=0
            elif 'name="ko:' in str(line) and 'entry' in str(line):
                Newentry = 1
                f1.write(str(line))
                if any(x in str(line).split('name="')[1].split('"')[0] for x in KO_meta):
                    Lable = 'ko'
                    try:
                        if any(x in str(line).split('reaction="')[1].split('"')[0] for x in Reac_meta):
                            Newtag = 1
                            try:
                                Reac_meta_filter.append(str(line).split('reaction="')[1].split('"')[0])
                            except IndexError:
                                pass
                    except IndexError:
                        Newtag = 1
            elif 'name="path:ko' in str(line) and 'entry' in str(line):
                if str(line).split(':ko')[1].split('"')[0] in Path_meta:
                    Newentry = 1
                    Lable = 'path'
                    f1.write(str(line))
            elif '<reaction' in str(line):
                Newentry=0
                Newcpd=0
                temp1 = str(line)
                temp2 = ""
                if any(x in str(line).split('name="')[1].split('"')[0] for x in Reac_meta):
                    Newcpd=1
                    Newentry = 1
            elif ('name="cpd:' in str(line) or 'name="gl:' in str(line)) and \
                    ('substrate' in str(line) or 'product' in str(line)):
                Lable = 'reaction'
                if str(line).split(':')[1].split('"')[0] not in Cpd_meta:
                    Newentry=0
                if Newcpd==1 and str(line).split(':')[1].split('"')[0] in Cpd_meta:
                    Cpd_meta_filter.append(str(line).split(':')[1].split('"')[0])
                temp2 = str(temp2) + str(line)
            elif '</reaction' in str(line) and Newentry==1:
                f1.write(temp1)
                f1.write(temp2)
                f1.write('    </reaction>\n')
                Newentry = 0
            elif "fgcolor" in str(line) and (Lable=='ko' or Lable=='cpd') and Newtag==1:
                f1.write(str(line).split('fgcolor="')[0] + 'fgcolor="#FF0000" ' + str(line).split('fgcolor="')[1].split(' ')[1])
            elif '</entry' in str(line):
                if Newentry ==1:
                    Newtag=0
                    Newentry=0
                    f1.write(str(line))
                else:
                    Newtag=0
                    Newentry = 0
            elif Newentry==1:
                f1.write(str(line))
        f1.write('</pathway>\n')
        f1.close()
        f1 = open(os.path.join(args.resultdir, str(os.path.splitext(str(File))[0])+'.ko01100.meta.filter.xml'), 'wb')
        Newentry = 1
        for line in open(os.path.join(args.resultdir, str(os.path.splitext(str(File))[0])+'.ko01100.meta.xml'), 'rb'):
            if ('name="cpd:' in str(line) or 'name="gl:' in str(line)) and 'entry' in str(line):
                Lable = 'cpd'
                if str(line).split(':')[1].split('"')[0] in Cpd_meta_filter:
                    f1.write(str(line))
                    Newentry = 1
                else:
                    Newentry = 0
            elif 'name="ko:' in str(line) and 'entry' in str(line):
                if any(x in str(line).split('name="')[1].split('"')[0] for x in KO_meta):
                    try:
                        if str(line).split('reaction="')[1].split('"')[0] in Reac_meta_filter:
                            Newentry = 1
                            f1.write(str(line))
                    except IndexError:
                        Newentry = 1
                        f1.write(str(line))
                    Lable = 'ko'
            elif 'name="path:ko' in str(line) and 'entry' in str(line):
                if str(line).split(':ko')[1].split('"')[0] in Path_meta:
                    Newentry = 1
                    Lable = 'path'
                    f1.write(str(line))
            elif '<reaction' in str(line):
                if any(x in str(line).split('name="')[1].split('"')[0] for x in Reac_meta_filter):
                    Newentry = 1
                else:
                    Newentry=0
                temp1 = str(line)
                temp2 = ""
            elif ('name="cpd:' in str(line) or 'name="gl:' in str(line)) and \
                    ('substrate' in str(line) or 'product' in str(line)):
                Lable = 'reaction'
                if str(line).split(':')[1].split('"')[0] not in Cpd_meta_filter:
                    Newentry = 0
                temp2 = str(temp2) + str(line)
            elif '</reaction' in str(line) and Newentry == 1:
                f1.write(temp1)
                f1.write(temp2)
                f1.write('    </reaction>\n')
                Newentry=0
            elif "fgcolor" in str(line) and (Lable == 'ko') and Newentry == 1:
                f1.write(
                    str(line).split('fgcolor="')[0] + 'fgcolor="#FF0000" ' + str(line).split('fgcolor="')[1].split(' ')[
                        1])
            elif '</entry' in str(line):
                if Newentry == 1:
                    f1.write(str(line))
                    Newentry=0
                else:
                    Newentry = 0
            elif Newentry == 1:
                f1.write(str(line))
        f1.write('</pathway>\n')
        f1.close()
    except IOError:
        print 'Files were missing?'

################################################### Programme #######################################################
for file_name in list_file:
    if Anno_type==0:
        Col_gene = int(args.gene) - 1
        Col_ko = int(args.ko) - 1
        Sept = str(args.format)
        Path=dict()
        for line in open(os.path.join(kegg_path, 'ko_map.txt'), 'rb'):
            Path.setdefault(str(line).split('\t')[0],str(line).split('\t')[1].split('\n')[0].split(' '))
        KO_list = dict()
        KO_list = intkolist(file_name, Col_gene, Col_ko, Sept)
        Path_list = dict()
        Path_list = pathassign(KO_list, Path)
        writepath(Path_list, KO_list)
        del file_name
    if Anno_type==1:
        pathlist(file_name)
