# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

std_table = ({
'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
'GGG': 'G' },
('TTG', 'CTG', 'ATG'),
('TAA', 'TAG', 'TGA')
)

oth_table= ({
'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
'GGG': 'G' },
('TTG', 'ATG'),
('TAA', 'TAG', 'TGA')
)            

def printCodonTable(table=std_table):
    new_table={}
    aminoSet=set(table[0].values())
    for s in aminoSet :
        new_table[s]=[]
    
    for item in table[0].items() :
        new_table[item[1]].append(item[0])
    
    print 'Tabulka obsahuje ',len(table[0].keys()),' trojic DNA, ktore reprezentuju ',len(aminoSet),' aminokyselin';
    print 'Aminokyselina',':',' Trojice baz';

    for item in new_table.items() :
        print '-----',item[0],'-----',':',item[1]
    
    print 'StartKodony su ',table[1];
    print 'EndKodony su ',table[2];
    
printCodonTable(oth_table)

# <codecell>

def cod2Aa(codon,table=std_table):
    upCodon=codon.upper()
    if table[2].count(upCodon)>0 :
        print codon,'je STOP kodon'
        return 'Stop';
    elif table[0].has_key(upCodon) :
        print 'Ku kodonu existuje aminokyselina ',table[0].get(codon);
        return table[0].get(codon)
    else : 
        print 'V tabulke sa takyto kodon nenachadza';
        return ''

# <codecell>

a=cod2Aa('TAA')

# <codecell>

def Aa2Cod(amino,table=std_table):
    if amino=='Stop' :
        print 'Stop kodony',table[2];
        return table[2]
    else :
        retCodons=[]
        for item in table[0].items() :
            if item[1]==amino :
                retCodons.append(item[0])
        print 'Pre aminokyselinu',amino,'existuju kodony',retCodons
        return retCodons

# <codecell>

a=Aa2Cod('F')

# <codecell>

def Dna2Cod(dna,pos=0,table=std_table):
    if (pos > len(dna)) | (pos < 0) :
        print 'Premenna pos je vacsia ako dlzka dna retazca';
        return
    
    dna=dna.upper()
    dnaRetazec=['~']
    for i in range(pos,len(dna),3) :
        aktualnyKodon=dna[i:i+3]
        if table[2].count(aktualnyKodon)>0 :   
            dnaRetazec.append('~')
            break
        elif table[0].has_key(aktualnyKodon) :
            dnaRetazec.append(table[0].get(aktualnyKodon))
        else : 
            dnaRetazec.append('_')
    return ''.join(dnaRetazec)

# <codecell>

Dna2Cod('TTDTTTTTTCAGCAACATGAGATTCTC',2)

# <codecell>


