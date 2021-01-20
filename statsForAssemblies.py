#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
import regex as re
import os,sys


def prepFolders():
    if os.path.isdir('temp')!=True:
        os.mkdir('temp')

def getFolders(inpath):
    files = os.path.listdir()

def loadGenome(infasta):
    return {i.description:str(i.seq) for i in SeqIO.parse(str(infasta),'fasta')}

def getLens(dct):
    return sorted([len(i) for i in dct.values()], reverse=True)

def calc_N50(lens):
    print('calculating N50 . . .')
    total = sum(lens)
    seq=0
    for i in lens:
        seq+=i
        if (seq/total)*1000000000 >= 0.5*1000000000:
            return i

def calc_NG50(lens, genome_estimate):
    print('calculating NG50 . . .')
    total = int(genome_estimate)
    seq=0
    for i in lens:
        seq+=i
        if (seq/total)*1000000000 >= 0.5*1000000000:
            return i

def max_minScaff(lens):
    return max(lens), min(lens)

def total_bp(scaffs):
    tot = 0
    for v in list(scaffs.values()):
        tot += len(v)
    return tot
    print('total size of assembly is:',tot/1000000,'mbp')

def checkGaps(assembly):
    Ns={k:[] for k in assembly.keys()}
    for k,v in assembly.items():
        for i in re.compile(r'N{5,}').finditer(v):
            Ns[k]+=[[i.start(), i.end()]]
    N={k:v for k,v in Ns.items() if len(v) >= 1}
    gaps=[]
    for k,v in N.items():
        for i in v:
            gaps.append(i[1]-i[0])
    return N, gaps

def telomerePrep():
    C4A2='CCCCAA'
    C3A3='CCCAAA'
    T2G4=str(Seq(C4A2).reverse_complement())
    T3G3=str(Seq(C3A3).reverse_complement())
    return C4A2, C3A3, T2G4, T3G3

def checkTelomeres(assembly):
    c4a2Telos=[]
    c3a3Telos=[]
    t2g4Telos=[]
    t3g3Telos=[]
    C4A2, C3A3, T2G4, T3G3 = telomerePrep()
    for k,v in assembly.items():
        if len(list(re.compile(r'('+C4A2+'){2,}').finditer(v[0:251]))) > 0:
            c4a2Telos.append(k)
        if len(list(re.compile(r'('+C3A3+'){2,}').finditer(v[0:251]))) > 0:
            c3a3Telos.append(k)
        if len(list(re.compile(r'('+T2G4+'){2,}').finditer(v[-251:-1]))) > 0:
            t2g4Telos.append(k)
        if len(list(re.compile(r'('+T3G3+'){2,}').finditer(v[-251:-1]))) > 0:
            t3g3Telos.append(k)
    c4a2Telos=list(set(c4a2Telos))
    c3a3Telos=list(set(c3a3Telos))
    t2g4Telos=list(set(t2g4Telos))
    t3g3Telos=list(set(t3g3Telos))
    fwdBoth=set(c4a2Telos+c3a3Telos)
    revBoth=set(t2g4Telos+t3g3Telos)
    bothTelos=list(fwdBoth & revBoth)
    singleTelos=list(fwdBoth ^ revBoth)
    return bothTelos, singleTelos

# def checkTelomeres(assembly):
#     c4a2Telos=[]
#     c3a3Telos=[]
#     t2g4Telos=[]
#     t3g3Telos=[]
#     C4A2, C3A3, T2G4, T3G3 = telomerePrep()
#     for k,v in assembly.items():
#         if len(list(re.compile(r'('+C4A2+'){2,}').finditer(v[251:-251]))) > 0:
#             c4a2Telos.append(k)
#         if len(list(re.compile(r'('+C3A3+'){2,}').finditer(v[251:-251]))) > 0:
#             c3a3Telos.append(k)
#         if len(list(re.compile(r'('+T2G4+'){2,}').finditer(v[251:-251]))) > 0:
#             t2g4Telos.append(k)
#         if len(list(re.compile(r'('+T3G3+'){2,}').finditer(v[251:-251]))) > 0:
#             t3g3Telos.append(k)
#     c4a2Telos=list(set(c4a2Telos))
#     print(len(c4a2Telos))
#     c3a3Telos=list(set(c3a3Telos))
#     print(len(c3a3Telos))
#     t2g4Telos=list(set(t2g4Telos))
#     print(len(t2g4Telos))
#     t3g3Telos=list(set(t3g3Telos))
#     print(len(t3g3Telos))
#     fwdBoth=set(c4a2Telos+c3a3Telos)
#     print(len(fwdBoth))
#     revBoth=set(t2g4Telos+t3g3Telos)
#     print(len(revBoth))
#     bothTelos=list(fwdBoth & revBoth)
#     print(len(bothTelos))
#     singleTelos=list(fwdBoth ^ revBoth)
#     print(len(singleTelos))
#     return bothTelos, singleTelos

def tblastn(inprot, assembly_db, evalue):
    out = 'temp/'+str(inprot).split('/')[-1]+'_v_'+str(assembly_db).split('/')[-1]
    os.system('tblasn -query '+str(inprot)+' -db '+str(assembly_db)+' -evalue '+str(evalue)+' -max_target_seqs 1 -outfmt 6 -out '+out)
    return out

def filterBlast(blrep, pattern):
    lt=[]
    for i in blrep:
        if i.startswith(str(pattern)):
            if i not in lt:
                lt.append(i)
    lst=list(set([i.split('\t')[0] for i in lt]))
    return lst

def getMappedPerc(blastRep,ref_fasta):
    CiliateProtDB = {i.description:str(i.seq) for i in SeqIO.parse(str(ref_fasta),'fasta')}
    cilProtPTET ={k:v for k,v in CiliateProtDB.items() if 'PTET' in k}
    cilProtPCAU={k:v for k,v in CiliateProtDB.items() if 'PCAU' in k}
    cilProtICH={k:v for k,v in CiliateProtDB.items() if 'IMG' in k}
    cilProtTTH={k:v for k,v in CiliateProtDB.items() if 'TTHERM' in k}
    blastRp = open(str(blastRep)).readlines()
    ptet=filterBlast(blastRp, 'PTET')
    pcau=filterBlast(blastRp, 'PCAU')
    ttherm=filterBlast(blastRp, 'TTHER')
    ich=filterBlast(blastRp, 'IMG')
    perc_ptet=len(ptet)/len(cilProtPTET)
    perc_pcau=len(pcau)/len(cilProtPCAU)
    perc_tthm=len(ttherm)/len(cilProtTTH)
    perc_ich=len(ich)/len(cilProtICH)
    return perc_ptet, perc_pcau, perc_tthm, perc_ich

if __name__ == "__main__":
    if len(sys.argv[1:]) != 6:
        print("You're missing some arguments!!\nUsage:\n\tpython statsForAssemblies.py [your assembly] []")
        sys.exit()

    assembly = sys.argv[1]

    evalue = sys.argv[]
