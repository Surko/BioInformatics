# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import subprocess as proc
import os
from itertools import combinations as comb

# <codecell>

def readFasta(fName, verbose=0) :
    
    startSeqCharacter = '>'
    retDictionary = {}
    sequence=''
    key=''
    try :
        fileObj = open(fName)
        for line in fileObj :            
            if len(line)==0 :                
                continue;     
            line = line.splitlines()[0]
            if line.startswith(startSeqCharacter) :                
                if key != '' :
                    retDictionary[key]=sequence
                sequence=''
                key=line.split(' ',1)[0][1:]                                  
                continue
            sequence+=line
        if key != '' :
            retDictionary[key]=sequence
    except IOError:
        print 'Can\'t open the file', fName
    finally :
        fileObj.close()
    print 'Done with verbose=',verbose    
    if verbose == 1 :
        for key, sequence in retDictionary.iteritems() :
            print 'The sequence',key,'with',len(sequence),'characters was read'
        
    return retDictionary

# <codecell>

# Vykreslenie Distance matice podla formatu
def printDistanceMatrix(title,keys,matrix) :
    # Dlzka Bunky v ktorej bude meno sekvencie
    length = max(map(len,keys)) + 2
    # Vytvaranie formatu {0:<x} - Prvy prvok odsadeny od lava s dlzkou bunky x
    string = ['{0:<',str(length),'}']
    # Vytvorenie formatu pre hodnoty v matici
    matrixformat = [ '{'+str(i)+':8.5f}' for i in range(1,len(matrix[0])+1)]    
    # Spojenie formatu
    string = string + matrixformat    
    # Vytvorenie textu s formatom
    string=''.join(string) 
    # Jeden riadok matice
    maineString ='';
    for i in range(0,len(keys)) :
        # Vytvorenie riadku pomocou funkcie format. Na vygenerovanie parametrov z listu pouzijeme *
        maineString=string.format(*([keys[i]]+matrix[i]))
        # V prvom riadku vypiseme meno kompresie a pocet sekvencii
        if i==0 :
            print (len(maineString)-len(title))/2 * ' ' + title
            print str(len(keys))
        print maineString
    # Nakonci vypiseme oddelovac matic
    print len(maineString)*'-'
    

#default directory
defDirectory = 'd:\\cvicenia\\BioInformatics\\du2\\';
#individual directories
indDirectory = 'd:\\cvicenia\\BioInformatics\\du2\\gencompress\\individual\\';
indzipDirectory = 'd:\\cvicenia\\BioInformatics\\du2\\zip\\individual\\'
#combinated directory
combDirectory = 'd:\\cvicenia\\BioInformatics\\du2\\gencompress\\combinated\\';
#concanated directories
concDirectory = 'd:\\cvicenia\\BioInformatics\\du2\\gencompress\\concanated\\';
conczipDirectory = 'd:\\cvicenia\\BioInformatics\\du2\\zip\\concanated\\';

# Precitanie suboru fasta s ulozenymi genomami
dic = readFasta(defDirectory + 'All_seq.fasta')
# Ziskanie klucov <=> mien sekvencii
keys = dic.keys()
# Usporiadanie mien podla abecedy seqA, seqB, ... atd. Nutne kvoli usporiadaniu do distance matrix. Indexy v tabulke priamo oznacuju v poradi co za sekvenciu to je. 
keys.sort()

# Dictionary pre dlzky skomprimovanych suborov podla kompresie GenCompress
individual = {}   
# Dictionary pre dlzky skomprimovanych suborov podla kompresie 7zip
individualzip = {}

# Ziskanie dlzok individualnych skomprimovanych suborov pomocou funkcie os.stat
for key in dic :        
    statInfo = os.stat(indDirectory + key + '.gen')
    individual[key]=statInfo.st_size
    statInfo = os.stat(indzipDirectory + key + '.7z')
    individualzip[key]=statInfo.st_size
    
# Dictionary pre dlzky skomprimovanych suborov spojenych sekvencii podla kompresie GenCompress
concanated = {}
# Dictionary pre dlzky skomprimovanych suborov spojenych sekvencii podla kompresie 7zip
concanatedzip = {}

# Ziskanie dlzok spojenych skomprimovanych suborov. For cyklus prechadza kombinaciami => v jednej iteracii mozme dosadit 2 dlzky.
for (fst,snd) in comb(keys,2) :
    statInfo = os.stat(concDirectory + fst + snd + '.gen')
    concanated[fst+snd]=statInfo.st_size
    statInfo = os.stat(concDirectory + snd + fst + '.gen')
    concanated[snd+fst]=statInfo.st_size
    statInfo = os.stat(conczipDirectory + fst + snd + '.7z')
    concanatedzip[fst + snd]=statInfo.st_size
    statInfo = os.stat(conczipDirectory + snd + fst + '.7z')
    concanatedzip[snd + fst]=statInfo.st_size

# Dictionary pre dlzky skomprimovanych suborov referencnych sekvencii podla kompresie GenCompress
combinated = {}
# Dictionary pre dlzky referencnych suborov spojenych sekvencii podla kompresie 7zip
combinatedzip = {}

# Ziskanie dlzok referencnych skomprimovanych suborov. For cyklus prechadza kombinaciami => v jednej iteracii mozme dosadit 2 dlzky.
for (fst,snd) in comb(keys,2) :
    statInfo = os.stat(combDirectory + fst + snd + '.gen')
    combinated[fst+snd]=statInfo.st_size 
    statInfo = os.stat(combDirectory + snd + fst + '.gen')
    combinated[snd+fst]=statInfo.st_size    
    # Pre zip subory musime pouzit vzorec pre aproximaciu K(u|v):=|Compress(vu)|-|Compress(v)|
    combinatedzip[fst+snd]=concanatedzip[snd+fst]-individualzip[snd]
    combinatedzip[snd+fst]=concanatedzip[fst+snd]-individualzip[fst]

# Distance Matica pre GenCompress kompresiu
distance = [len(keys) * [0.0] for i in range(len(keys))]
# Distance Matica pre 7zip komrpesiu
distancezip = [len(keys) * [0.0] for i in range(len(keys))]

# Vypisanie Distance Matic podla aproximacie Kolmogorov Complexity D(u,v) = 1 - R(u,v), kde R(u,v) je podla slajdov. 
for i in range(len(keys)) :
    for j in range(len(keys)) :
        if i == j:
            continue
        else :            
            # distance pracuje s GenCompress
            distance[i][j]=round(1.0-(float)((individual[keys[i]] - combinated[keys[i]+keys[j]]))/concanated[keys[i]+keys[j]],5)            
            # distance pracuje so 7zip
            distancezip[i][j] = round(1.0-(float)((individualzip[keys[i]] - combinatedzip[keys[i]+keys[j]]))/concanatedzip[keys[i]+keys[j]],5)            

# Prvy sposob vykreslenia matic --- neprehladny            
#print 'Distance with GenCompress' 
#print '---------------------'
#print distance
#print 'Distance with 7z' 
#print '---------------------'
#print distancezip

# Vykreslenie Distance Matic pre jednotlive kompresie. Vykreslenie je vykonane funkciou printDistanceMatrix, ktora vykresli matice vo formate, ktory je podporovany programami na generovanie Newikovych a fylogenetickych stromov.
printDistanceMatrix('GenCompress',keys,distance)
printDistanceMatrix('7z',keys,distancezip)

# Vypisanie info o strankach kde sa daju vygenerovat Newickove stromy a Fylogeneticke stromy
print ' Teraz mozme pouzit vytvorenie Newickovho stromu napriklad zo stranky http://mobyle.pasteur.fr/cgi-bin/portal.py?#forms::bionj'
print ' Po vygenerovani newickovho stromu mozme zase vyuzit stranku http://iubio.bio.indiana.edu/treeapp/treeprint-form.html na vygenerovanie fylogenetickeho stromu'

# <codecell>

# Vytvorenie suborov so sekvenciami ako individualnych tak aj spojenych

# Zakladna zlozka kde sa nachadza GenCompress program aj so subormi na kompresiu
defDirectory = 'd:\\cvicenia\\BioInformatics\\du2\\';

# Precitanie suboru fasta s ulozenymi genomami
dic = readFasta(defDirectory + 'All_seq.fasta')
# Ziskanie klucov <=> mien sekvencii
keys = dic.keys()
# Usporiadanie mien podla abecedy seqA, seqB, ... atd. Neni nutne
keys.sort()

# Vytvorenie individualnych sekvencii. Vyberie sekvenciu z dictionary, v ktorej su ulozene vsetky sekvencie a zapise do suboru s menom popisujuci sekvenciu.
for key in keys :
    myFile = open(defDirectory + key,'w');
    myFile.write(dic[key])
    myFile.close()

# Vytvorenie spojenych sekvencii. Vyberie prvu aj druhu sekvenciu z dictionary, v ktorej su ulozene vsetky sekvencie a zapise ich za seba do suboru s menom popisujuci spojenu sekvenciu.
for (fst,snd) in comb(keys,2) :
    myFile = open(defDirectory + fst + snd,'w');
    myFile.write(dic[fst] + dic[snd])
    myFile.close()
    myFile = open(defDirectory + snd + fst,'w');
    myFile.write(dic[snd] + dic[fst])
    myFile.close()

# <codecell>

# Vytvorenie individualnych suborov (A)

# Definovane premenne so zlozkami v ktorych sa nachadzaju subory na skomprimovanie

# Zakladna zlozka kde sa nachadza GenCompress program aj so subormi na kompresiu
defDirectory = 'd:\\cvicenia\\BioInformatics\\du2\\';
# Zlozka kde sa nachadza 7z program.
zipFolder = 'd:\\utils\\Totalcmd\\arch\\7-Zip\\';
# Zlozka kam sa ulozia skomprimovane subory pomocou kompresie GenCompress
ind1Directory = 'd:\\cvicenia\\BioInformatics\\du2\\gencompress\\individual\\';
# Zlozka kam sa ulozia skomprimovane subory pomocou kompresie 7z
ind2Directory = 'd:\\cvicenia\\BioInformatics\\du2\\zip\\individual\\';
# Precitanie suboru fasta s ulozenymi genomami
dic = readFasta(defDirectory + 'All_seq.fasta')
# Ziskanie klucov <=> mien sekvencii
keys = dic.keys()

for key in dic :  
    # Komprimovanie suborov pomocou GenCompress do zlozky defDirectory
    proc.call([defDirectory + 'GenCompress.exe', defDirectory + key]); 
    # Presunutie vygenerovane suboru do zlozky pre individualne subory
    os.renames(defDirectory + key + '.gen',ind1Directory + key + '.gen')   
    # Komprimovanie suborov pomocou 7z do zlozky defDirectory
    proc.call([zipFolder + '7z.exe', 'a', defDirectory + key + '.7z', defDirectory + key]);
    # Presunutie vygenerovanych suborov do zlozky pre individualne subory
    os.renames(defDirectory + key + '.7z',ind2Directory + key + '.7z');
    # Upratanie
    os.remove(defDirectory + key + '.log')

# <codecell>

# Vytvorenie kombinovanych suborov (pomocou referencie A|R)

# Definovane premenne so zlozkami v ktorych sa nachadzaju subory na skomprimovanie

# Zakladna zlozka kde sa nachadza GenCompress program aj so subormi na kompresiu
defDirectory = 'd:\\cvicenia\\BioInformatics\\du2\\';
# Zlozka kde sa nachadza 7z program.
zipFolder = 'd:\\utils\\Totalcmd\\arch\\7-Zip\\';
# Zlozka kam sa ulozia skomprimovane subory pomocou kompresie GenCompress
comb1Directory = 'd:\\cvicenia\\BioInformatics\\du2\\gencompress\\combinated\\';
# Zlozka kam sa ulozia skomprimovane subory pomocou kompresie 7z. Nepouzite lebo pri 7z nevytvarame refenciou
comb2Directory = 'd:\\cvicenia\\BioInformatics\\du2\\zip\\combinated\\';
# Precitanie suboru fasta s ulozenymi genomami
dic = readFasta(defDirectory + 'All_seq.fasta')
# Ziskanie klucov <=> mien sekvencii
keys = dic.keys()
# Usporiadanie mien podla abecedy seqA, seqB, ... atd. Neni nutne
keys.sort()  

for fst in keys :
    for snd in keys :
        # Ked su subory rovnake tak nekomprimujem
        if fst == snd :
            continue
        # Komprimovanie referenciou na iny subory pomocou GenCompress do zlozky defDirectory
        proc.call([defDirectory + 'GenCompress.exe', defDirectory + fst,'-c',defDirectory + snd])
        # Presunutie vygenerovane suboru do zlozky pre kombinovane subory
        os.renames(defDirectory + fst + '.gen',comb1Directory + fst + snd + '.gen')    
        # Upratanie
        os.remove(defDirectory + fst + '.log')

# <codecell>

# Vytvorenie spojenych suborov (AB)

# Definovane premenne so zlozkami v ktorych sa nachadzaju subory na skomprimovanie

# Zakladna zlozka kde sa nachadza GenCompress program aj so subormi na kompresiu
defDirectory = 'd:\\cvicenia\\BioInformatics\\du2\\';
# Zlozka kde sa nachadza 7z program.
zipFolder = 'd:\\utils\\Totalcmd\\arch\\7-Zip\\';
# Zlozka kam sa ulozia skomprimovane subory pomocou kompresie GenCompress
conc1Directory = 'd:\\cvicenia\\BioInformatics\\du2\\gencompress\\concanated\\';
# Zlozka kam sa ulozia skomprimovane subory pomocou kompresie 7z.
conc2Directory = 'd:\\cvicenia\\BioInformatics\\du2\\zip\\concanated\\';
# Precitanie suboru fasta s ulozenymi genomami
dic = readFasta(defDirectory + 'All_seq.fasta')
# Ziskanie klucov <=> mien sekvencii
keys = dic.keys()
# Usporiadanie mien podla abecedy seqA, seqB, ... atd. Neni nutne
keys.sort()  

for fst in keys :
    for snd in keys :
        # Ked su subory rovnake tak nekomprimujem
        if fst == snd :
            continue  
        # Komprimovanie suborov pomocou GenCompress do zlozky defDirectory
        proc.call([defDirectory + 'GenCompress.exe', defDirectory + fst + snd])
        # Presunutie vygenerovanych suborov do zlozky pre spojene subory
        os.renames(defDirectory + fst + snd + '.gen',conc1Directory + fst + snd + '.gen')
        # Komprimovanie suborov pomocou 7z do zlozky defDirectory
        proc.call([zipFolder + '7z.exe', 'a', defDirectory + fst + snd + '.7z', defDirectory + fst + snd]);
        # Presunutie vygenerovanych suborov do zlozky pre spojene subory
        os.renames(defDirectory + fst + snd + '.7z',conc2Directory + fst + snd + '.7z');
        # Upratanie
        os.remove(defDirectory + fst + snd + '.log')

