import os
import argparse
import copy
import math

############################################ Arguments and declarations ##############################################
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i",
                    help="A single file name", type=str, default='PAO_KO.txt',
                    metavar='your file with KO numbers and their color lable (integer)')
parser.add_argument("-ko",
                    help="The colomn number of KO", type=int, default=1,
                    metavar=1)
parser.add_argument('-col',
                    default="2", action='store', type=str,
                    metavar='2,3,4',
                    help="The colomns number of expression levels (separated by ',')")
parser.add_argument("--f",
                    help="The seperation character of input file (default: tab)", type=str, default="\t",metavar=
                    "tab: default or $'\t' , blackspace: ' '")
parser.add_argument("--p",
                    help="The method for merging labels for path with multiple subpaths "
                         "(0 for common, 1 for max, 2 for min)",
                    type=int, default=2,metavar='0 or 1 or 2',choices=[0, 1, 2])
parser.add_argument("--v",
                    help="CRPKM ranks or CRPKM value (default)"
                         "(r for CRPKM ranks, v for CRPKM value)",
                    type=str, default='v',metavar='r or v',choices=['r', 'v'])
################################################## Definition ########################################################
args = parser.parse_args()
input_path = os.path.abspath(args.i)
file_pao = os.path.join(input_path)
filedir, filename = os.path.split(file_pao)
KO_col=int(args.ko)-1
Lable_col=[]
for col in args.col.split(','):
    Lable_col.append(int(col.strip())-1)
sept = str(args.f)
try:
    os.mkdir('output')
except OSError:
    pass

################################################### Function #######################################################
def setpath(Number,EC,list1,list2):
    if EC=='None':
        list1[Number]=[['None']]
    elif EC=='No':
        list1[Number]=[['No']]
    elif 'and' in EC:
        Position=0
        list1.setdefault(Number,[])
        for key in EC.split('and'):
            setko(Number,key,list2,Position)
            Position += 1
            list1[Number].append([])
    else:
        list1.setdefault(Number, [[]])
        setko(Number, EC, list2, 0)


def setko(Number,EC,list2,Position):
    if 'or' in EC:
        for key in EC.split('or'):
            if list2.get(key.strip(),'None')=='None':
                list2.setdefault(key.strip(),[str(Number) + ':' + str(Position)])
            elif Number!=list2.get(key.strip(),'None')[0]:
                list2[key.strip()].append(str(Number) + ':' + str(Position))
    else:
        if EC.strip()!='':
            if list2.get(EC.strip(), 'None') == 'None':
                list2.setdefault(EC.strip(), [str(Number) + ':' + str(Position)])
            elif Number != list2.get(EC.strip(), 'None')[0]:
                list2[EC.strip()].append(str(Number) + ':' + str(Position))


def tablemapping(KO,list1,list2,lable,level):
    KOnew=''
    if lable >0:
        lable=int((math.log10(lable)-math.log10(minvalue))/(math.log10(maxvalue)-math.log10(minvalue))*level)
    if KO in list2:
        KOnew=KO.capitalize()
    elif KO.capitalize() in Listko:
        KOnew=Listko[KO.capitalize()]
    if KOnew !='':
        if KOnew in list2:
            for key in list2[KOnew]:
                list1[key.split(':')[0]][int(key.split(':')[1])].append(lable)




def pathmerging(list1,method,f):
    for number in list1:
        if len(list1[number])>1:
            position = 0
            for key in list1[number]:
                list1[number][position] = comparison(key, 1)  # for isoenzymes
            list1[number] = comparison(list1[number], method)
        else:
            list1[number] = comparison(list1[number][0], 1)  # for isoenzymes
        f.write(str(number)+'\t'+str(list1[number])+'\n')


def comparison(lablelist,method):
    if len(lablelist)==0:
        return 0
    elif len(lablelist)==1:
        return lablelist[0]
    elif len(lablelist)>1:
        if method==1:
            return max(lablelist)
        elif method==2:
            return min(lablelist)
        elif method==0:
            frequency={x:lablelist.count(x) for x in lablelist}
            lables, fres = frequency.keys(), frequency.values()
            if len(fres.index(max(fres)))==1:
                return lables[fres.index(max(fres))]
            else:
                return sum(lables[key] for key in fres.index(max(fres)))/len(fres.index(max(fres)))

################################################### Programme #######################################################
#initiation
Listko=dict()
for line in open('KO_EC_mapping.txt','rb'):
    if str(line).split('\t')[3].split('\r')[0].split('\n')[0]!='':
        try:
            Listko.setdefault(str(line).split('\t')[1].split('ko:')[1].strip().capitalize(),
                          str(line).split('\t')[3].split('ec:')[1].split('\r')[0].split('\n')[0].strip())
        except IndexError:
            print 'KO_EC: '+str(line).split('\t')[1].split('ko:')[1].strip().capitalize()


Listnode=dict()
KOnode=dict()
for line in open('Node_number.txt','rb'):
    if str(line).split('\t')[1].split('\r')[0].split('\n')[0]!='':
        setpath(str(line).split('\t')[0],str(line).split('\t')[1].split('\r')[0].split('\n')[0],Listnode,KOnode)
Listnode1=copy.deepcopy(Listnode)
KOnode1=copy.deepcopy(KOnode)


Listedge=dict()
KOedge=dict()
for line in open('Edge_number.txt','rb'):
    if str(line).split('\t')[1].split('\r')[0].split('\n')[0]!='':
        setpath(str(line).split('\t')[0],str(line).split('\t')[1].split('\r')[0].split('\n')[0],Listedge,KOedge)
Listedge1=copy.deepcopy(Listedge)
KOedge1=copy.deepcopy(KOedge)

#set the maxmum and minimum expression value
maxvalue=minvalue=0.0
for line in open(filename, 'rb'):
    KO = str(line).split(sept)[KO_col].strip()
    if KO in KOnode or KO in KOedge or Listko.get(KO,'None') in KOnode \
            or Listko.get(KO,'None') in KOedge:
        for col in Lable_col:
            try:
                if float(str(line).split(sept)[col].strip()) >0:
                    maxvalue=max(maxvalue,float(str(line).split(sept)[col].strip()))
                    if minvalue ==0.0:
                        minvalue=float(str(line).split(sept)[col].strip())
                    else:
                        minvalue = min(minvalue, float(str(line).split(sept)[col].strip()))
            except ValueError:
                pass

#mapping expresison value to node and edge tables by ko number
if maxvalue==0:
    print 'No valid expression value (should be >0)!'
else:
    for col in Lable_col:
        Listnode = copy.deepcopy(Listnode1) #reset node and edge table
        KOnode = copy.deepcopy(KOnode1)
        Listedge = copy.deepcopy(Listedge1)
        KOedge = copy.deepcopy(KOedge1)
        f1 = open('output/{0}_{1}_edge.txt'.format(str(os.path.splitext(filename)[0]), str(col)), 'ab')
        f2 = open('output/{0}_{1}_node.txt'.format(str(os.path.splitext(filename)[0]), str(col)), 'ab')
        f1.write('No\tColor\n')
        f2.write('No\tColor\n')
        for line in open(filename, 'rb'):
            KO = str(line).split(sept)[KO_col].strip()
            if args.v=='r':
                try:
                    tablemapping(KO, Listnode, KOnode, float(str(line).split(sept)[col].strip()),3)
                    tablemapping(KO, Listedge, KOedge, float(str(line).split(sept)[col].strip()),3)
                except ValueError:
                    pass
            elif args.v == 'v':
                try:
                    tablemapping(KO, Listnode, KOnode, float(str(line).split(sept)[col].strip()),9)
                    tablemapping(KO, Listedge, KOedge, float(str(line).split(sept)[col].strip()),9)
                except ValueError:
                    pass
            else:
                print 'Invalid input of --v parameter: please input r for CRPKM ranks, v for LCRPKM value.'
        pathmerging(Listedge, args.p, f1)
        pathmerging(Listnode, args.p, f2)
        f1.close()
        f2.close()


