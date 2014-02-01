# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Gibbs Sampling

# <markdowncell>

# V súbore MotSam.fasta máme 10 DNA sekvencií dĺžky 600 báz, v ktorých je (umelo pridaný) motív dĺžky 15 symbolov s maximálne 4 mutáciami v každom reťazci. Naprogramujte vlasntú verziu Gibbsovho vzorkovania a pokúste sa nájsť tento motív. Odovzdajte
# Program - vašu implementáciu Gibbsovho vzorkovania. Navrhnite, čo by mal program zobrazovať v priebehu výpočtu, aby sme mohli sledovať vývoj kvality nájdených motívov.
# Textový popis vášho riešenia s výsledkami 
# aký motív ste našli s akými parametrami programu v sekvenciách v súbore MotSam.fasta,
# ako dlho to trvalo, čo bolo vidieť v priebehu výpočtu.
# Pre nájdený motív uveďte jeho skóre, tj. vzdialenosť k reťazcom DNA.
# Ako by sa dalo vylepšiť Gibbsovo vzorkovanie, aby bolo spoľahlivejšie, rýchlejšie a pod.
# Aké výsledky dáva váš program pri hľadaní motívu dĺžky 9 báz s maximálne 4 mutáciami na sekvenciách v súbore MotSamShort.fasta? Ako sa líšil priebeh výpočtu nad druhým zadaním od výpočtu nad prvým zadaním?

# <headingcell level=2>

# Pomocné funkcie pre gibbsovo vzorkovanie

# <codecell>

import random as r
import numpy as np
import os
import time

seq1file='d:\cvicenia\BioInformatics\du3\MotSam.fasta'
seq2file='d:\cvicenia\BioInformatics\du3\MotSamShort.fasta'

# <codecell>

def readFasta(fName, verbose=1) :
    
    startSeqCharacter = '>'
    retDictionary = {}
    sequence=''
    key=''
    try :
        fileObj = open(fName)
        for line in fileObj :            
            if len(line)==0 :                
                continue;
            if line.startswith(startSeqCharacter) :                
                if key != '' :                    
                    retDictionary[key]=sequence.partition('\n')[0]
                sequence=''
                key=line.split(' ',1)[0][1:]  
                continue
            sequence+=line
        if key != '' :            
            retDictionary[key]=sequence.partition('\n')[0]
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

# definovanie indexov pre baze
nuclMatrixInd={'a':0,'A':0,'c':1,'C':1,'g':2,'G':2,'t':3,'T':3};
nuclMatrixInv=['A','C','G','T']

def gibbsSampling(sequences,l,mut,N,choice=0,verbose=0) :
    size = len(sequences)
    
    # Ziskanie nahodnych l-merov zo sekvencii. Z kazdej sekvencie jeden.
    lMers=randInitMers(sequences,l,choice)    
        
    # Ziskanie 4 prvkoveho pola kde je vycet baz zo sekvencii.[#A,#C,#G,#T]
    allQ = nuclCount(sequences)
    
    # aktualna iteracia
    iteration = 0;
    # aktualne skore
    score = 0        
    
    #Predchadzjuce lmery - vyuzite pri moznosti 2
    prevlMers=np.empty((size,l),dtype=str)
    
    # Nastavenie najlepsieho motifu a najlepsieho skore pomocou funkcie getScore
    bestMotif,scoreBest = getScore(lMers,l)           
    # Nastavenie najlepsich l-merov do premennej bestMers
    bestMers = lMers.copy()
        
    # Beh algoritmu dokym neprebehne pocet iteracii N
    for i in range(N) :                                        
        if choice == 2 :
            if (prevlMers==lMers).all() :
                bestMers = lMers;
                bestMotif = motif;
                scoreBest = scoreOrig;
                break
            else :
                prevlMers = lMers.copy()
        
        # Nahodne cislo 
        h = r.randint(0,size-1);
        
        # Vystavenie profilu z a1,a2,...,ah-1,ah+1,...,at vybranych l-merov
        X = buildProfile(lMers[range(h)+range(h+1,size),:],l);
                
        # Nahodne vybrata sekvencia 
        DNAh = sequences[h];
        
        # Ziskanie poctu baz dokopy v sekvenciach
        Q = (allQ - nuclCount(DNAh));
        # A znich relativne cetnosti
        Q /= sum(Q)
        
        # Pravdepodobnostne pole w pomocou list comprehension
        w=[probW(DNAh[i:i+l],X,Q) for i in range(len(DNAh)-l)]
        # suma vsetkych w
        sumW = sum(w)
        
        # Ked je suma rovna nule tak pravdepodobnost vsetkych l-merov bola nulova.
        # V takomto pripade nastavime nove ah na nahodny l-mer a dalej iterujeme.
        if sumW == 0 :
            newAhIndex = r.randint(0,len(DNAh)-l)
        else :
            # Musime ziskat urcite a. Kazde a ma pravdepodobnost urcenu z w(a)/sum(w)
            w/=sumW;
            # Nahodne cislo na ziskanie a-cka. Kazda pravdepodobnost z w obsadzuje urcity vysek a podla pravdepodobnosti urcime, co za usek to je.            
            rProb = r.random();
            p = 0;
            newAhIndex = 0;
            
            while True :
                p += w[newAhIndex];
                if p > rProb :
                    break;
                newAhIndex+=1;                            
        
        # l-mer ah nastavime na novy. Z DNAh vyberieme sekvenciu dlzky l
        lMers[h,:]=list(DNAh[newAhIndex:newAhIndex+l])

        # Iterovanie
        iteration+=1
        
        #Motif a skore z aktualnych lmerov.
        motif,scoreOrig = getScore(lMers,l)
        
        # Ked je skore novych l-merov menej ako scoreBest tak su to zatial najlepsie lmery
        if scoreOrig < scoreBest :
            bestMotif = motif
            bestMers = lMers.copy()
            scoreBest = scoreOrig
           
    if verbose==1 :    
        print 'Vypocet som skoncil na',iteration,'iteracii' 
        
    return bestMotif,scoreBest,bestMers       
        

# <codecell>

def randInitMers(sequences,l,choice=0) :    
    size = len(sequences)
    # l-mery budeme mat ulozene v numpy matici stringov.
    lMers=np.empty((size,l),dtype=str)
    for i in range(size) :
        while True :     
            # Vybranie nahodnych l-merov z kazdej sekvencie
            newSeq = randomSelectLMer(sequences[i],l)            
            if choice==0 :
                # Pridavanie l-meru do riadku. Nato aby v kazdej bunke bol iba jeden char tak musim pouzit funkciu list
                lMers[i,:]=list(newSeq);
                break;
            elif choice==1 :
                # Just a test
                # Vyberanie iba takych sekvencii, ktore maju od ostatnych vzdialenost najviac 2 * mut    
                # Dangerously slow!!!
                dist=[hamDist(newSeq,_seq) for _seq in lMers]                                              
                if dist==[] or max(dist) <= 2 * mut :
                    lMers[i,:]=list(newSeq);
                    break; 
            else :
                lMers[i,:]=list(newSeq);
                break; 
    return lMers

# <codecell>

# Ziskanie skore a consensus motivu pre l-mery dlzky l
def getScore(lMers,l) :        
    score = 0
    # Pole kam budeme ukladat pismenka motivu
    bestMotif=[]
    for i in range(lMers.shape[1]) :
        # Pole ktore reprezentuje jeden stlpec skorovacej tabulky => 1. prvok je pocet A-cok v sekvenciach na 1. mieste, 2. prvok je pocet C-cok... atd        
        dictionar=[0,0,0,0]
        for c in lMers[:,i] :
            dictionar[nuclMatrixInd[c]]+=1
        # Po iteraciach mozme zistit ktora baza sa tam nachadza najviac => priradenie bazy do bestMotif. Skore ziskame sum(dictionar)-max(dictionar)
        maximum=max(dictionar)
        # Index maxima na zistenie pismenka
        maxIndex=dictionar.index(maximum)
        # Scitavanie skore pre kazdy stlpec
        score+=sum(dictionar)-maximum
        # Pridanie pismenka do bestMotif
        bestMotif.append(nuclMatrixInv[maxIndex])
    return ''.join(bestMotif),score

# <codecell>

def hamDist(seqA,seqB) :   
    #print seqA,seqB
    dist = 0;  
    for i in range(len(seqA)) :
        if seqA[i]!=seqB[i] :
            dist+=1
    return dist;

# <codecell>

# Ziskanie w(a) = P(a|X) / P(a|Q)
def probW(a,X,Q) :
    pX = 1.0;
    for i in range(len(a)) :
        pX*=X[nuclMatrixInd[a[i]],i]
    pQ = 1.0;
    for i in range(len(a)) :
        pQ*=Q[nuclMatrixInd[a[i]]]
    return pX/pQ
    

# <codecell>

# Metoda, ktora navrati nahodne vybrany lmer zo sekvencie seq
def randomSelectLMer(seq,l,verbose=0) :
    # Startovacia pozicia nahodne zvoleneho l-meru
    lmerStart = r.randint(0,len(seq)-l)
    if verbose==1 :
        print lmerStart
    # Ziskanie podretazca
    return seq[lmerStart:lmerStart+l]    

# <codecell>

# Pocitanie globalnych frekvencii bazi
def nuclCount(sequences) :
    nuclArr=np.array([0.0,0.0,0.0,0.0]);
    # Pre kazdu sekvenciu
    for seq in sequences :
        # Prechadzam kazdym pismenom a pridavam do dictionary
        for c in seq :
            nuclArr[nuclMatrixInd[c]]+=1
    return nuclArr

# <codecell>

# Vytvorenie profilu zo sekvencii zadanymi v poli sequencii. Dlzka profilu bude l
def buildProfile(sequences,l,verbose=0) :
    numOfSeq = len(sequences) + 4;
    # Preddefinovanie profilu, tak ze defaultne hodnoty budu 1. To znamena ze kazda baza sa nachadza ako keby o jednu viac
    profile = np.ones((4,l),dtype=float);
    # Vytvaranie poctoveho profilu pre kazdu sekvenciu. profile[baza,i] = Pocet bazi na i-tom mieste v danych sekvenciach
    for seq in sequences :
        for i in range(l) :
            profile[nuclMatrixInd[seq[i]],i]+=1;
    
    # Vytvorenie profilu relativneho poctu. profile[baza,i] = Pravdepodobnost ze na i-tom mieste v danych sekvenciach bola dana baza
    for i in range(l) :
        # Ziskanie relativneho poctu pre kazdu bazu => predelenie kazdej bunky v jednom stlpci sumou.
        profile[:,i]/=numOfSeq
    
    return profile

# <headingcell level=2>

# Hlavné vykonanie gibbsovho vzorkovania pre MotSam

# <markdowncell>

# Pustanie gibbsovho vzorkovania na MotSam.fasta s poctom baz 15 a maximalne 4 mutaciami. Do premennych motif,score,array ukladame navratove hodnoty z metody.Tie zodpovedaju consensus motivu, skore consensus motivu od ostatnych sekvencii => sucet vzdialenosti od ostatnych sekvencii a vypis najdenych sekvencii. 
# 
# Kedze skore obsahuje sucet vzdialenosti consensus motivu od ostatnych sekvencii a vieme ze sa maju vyskytovat najviac 4 mutacie, tak jedini vhodni kandidati su taki, kde skore bude menej ako 40 (pocetsekvencii * pocetmutacii). Nie kazdy takyto je ale optimalny. Tych najlepsich kandidatov si ukladame do 10 prvkoveho pola aby sme mohli vidiet ako sa vyvijalo hladanie riesenia.
# 
# Celkova rychlost gibbsovho vzorkovania zavisi od poctu bazi (teraz 15) ako aj od poctu iteracii v ramci jedneho vzorkovania(v mojom pripade 200) ci pocet opakovani gibbsovho vzorkovania (v mojom pripade 50).
# 
# Taktiez som skusal zvacsit pocet opakovani gibbsovho vzorkovania na 500.

# <codecell>

# Nacitanie sekvencii z fasta suboru
sequences1 = readFasta(seq1file,verbose=0);
# ziskanie hodnot
values = sequences1.values()
# a klucov, ktore nevyuzijeme
keys = sequences1.keys()

# Vytvorenie pola do ktoreho budeme ukladat 10 najlepsich behov
best=[['',150,[]] for i in range(10)]

timeArray=[]

# 50 krat opakujeme gibbsovo samplovanie s 
for j in range(50) :
    # Nastavime casovac pre kazdu iteraciu
    start=time.clock()
    # Samplovanie pustame s dlzkou motivu 15. Argument 200 znaci kolko krat opakujeme iterovanie motivmi
    motif,score,array = gibbsSampling(values,15,4,200,choice=0)
    # Vypisanie aktualnej iteracie opakovania a ziskaneho skore
    print 'Iter:',j,'with score',score,'and motif',motif
    # Nasledne priradim tuto trojicu do pola best na index kam patri. Ked pole pretecie tak vymazem posledny prvok (ten najvacsi)
    for k in range(10) :
        if score < best[k][1] :
            best.insert(k,[motif,score,array])
            del best[-1]
            break
    timeArray.append(time.clock()-start)        

# Vypis ako rychlo zbehol vypocet
print 'Gibbsovo vzorkovanie zbehlo v case',sum(timeArray)    

# <codecell>

# Nacitanie sekvencii z fasta suboru
sequences1 = readFasta(seq1file,verbose=0);
# ziskanie hodnot
values = sequences1.values()
# a klucov, ktore nevyuzijeme
keys = sequences1.keys()

# Vytvorenie pola do ktoreho budeme ukladat 10 najlepsich behov
best=[['',150,[]] for i in range(10)]

timeArray=[]

# 50 krat opakujeme gibbsovo samplovanie s 

for j in range(50) :
    # Nastavime casovac pre kazdu iteraciu
    start=time.clock()
    # Samplovanie pustame s dlzkou motivu 15.Choice 2 znamena ze algoritmus bezi dokym nebudu motivy v poli rovnake po kazdej iteracii
    motif,score,array = gibbsSampling(values,15,4,200,choice=2)
    # Vypisanie aktualnej iteracie opakovania a ziskaneho skore
    print 'Iter:',j,'with score',score,'and motif',motif
    # Nasledne priradim tuto trojicu do pola best na index kam patri. Ked pole pretecie tak vymazem posledny prvok (ten najvacsi)
    for k in range(10) :
        if score < best[k][1] :
            best.insert(k,[motif,score,array])
            del best[-1]
            break
    timeArray.append(time.clock()-start)        

print 'Gibbsovo vzorkovanie zbehlo v case',sum(timeArray)           

# <headingcell level=3>

# Vysledky na MatSam

# <codecell>

# Tento skript vyberie kazdu trojicu v ktorej su ulozene najlepsie najdene motify.
# Z kazdej trojice vyberie pole (na druhom mieste) v ktorych su kusky sekvencii a 
# pre kazde zistime hamming vzdialenost od consensus motifu. Ked su vzdialenosti od consensus motifu
# mensie alebo rovne hodnote mut tak je to spravny motif

mut = 4

for i in range(len(best)) :
    w=[hamDist(best[i][0],seq) for seq in best[i][2]]    
    if all([x <= mut for x in w]) :
        print 'Optimal',w,'Consensus motif',best[i][0]
    else :
        print w,'Consensus motif',best[i][0]

# <markdowncell>

# Vyššie je vidieť 10 najlepších kandidátov. Pri kazdom behu vyzera toto pole inak. Gibbsovo vzorkovanie nemusi vzdy nast akurat optimalne riesenie. Napr. teraz ma optimalne skore=25, niekedy akurat tento pripad nemusi najst, alebo v niektorych pripadoch moze najst tiez riesenie so skore=25, ktore by malo vyhovovat, no niektora sekvencia bude mat viac nez 4 mutacie.
# 
# V tomto pripade je prvý,druhy aj treti ten optimálny  a ako je vidieť vzdialenosť najdeneho motivu [GTGACACATTATGTC] od najdenych sekvencii je [3, 4, 3, 2, 2, 4, 0, 4, 3, 0] => najviac 4 mutácie. V ďalšom prípade je medzi jedným motivom a sekvenciou počet mutácii rovný 5.
# 
# Pri 500 behoch Gibbsovho vzorkovania sa nam matica s najlepšími kandidatmi vytriedila iba k jednemu consensus motivu, ktorý ma vzdialenosti od nájdených 10 motívov rovné [3, 4, 3, 2, 2, 4, 0, 4, 3, 0], co je rovnake vlastne riesenie ake sme nasli pri mensich poctoch opakovani

# <codecell>

best

# <markdowncell>

# Vypis obsahu najlepších 10 nájdených kanditátov. Každý kandidát je trojprvkové pole, kde prvy prvok je consensus motiv => najbližší motiv k sekvenciám uložených v treťom prvku. V druhom prvku je skóre => sčítané hammingovské vzdialenosti consensus motívu od sekvencií. 

# <markdowncell>

# \begin{align}
# VysledneSkore &= 25 \\
# VyslednyMotif &= \begin{array}{ccccccccccccccc} G & T & G & A & C & A & C & A & T & T & A & T & G & T & C \end{array} \\
# VysledneSekvencie &= \begin{pmatrix}
# G & A & G & A & C & T & C & A & T & A & A & T & G & T & C \\
# G & T & T & A & G & A & C & A & A & T & A & T & G & T & T \\
# C & T & G & A & A & A & C & A & C & T & A & T & G & T & C \\
# C & T & G & A & C & A & C & G & T & T & A & T & G & T & C \\
# G & G & G & A & C & A & C & A & T & A & A & T & G & T & C \\
# G & T & G & A & G & A & G & A & C & G & A & T & G & T & C \\
# G & T & G & A & C & A & C & A & T & T & A & T & G & T & C \\
# G & T & G & A & T & A & T & G & T & T & A & T & G & A & C \\ 
# G & T & G & A & C & A & C & A & T & G & G & G & G & T & C \\
# G & T & G & A & C & A & C & A & T & T & A & T & G & T & C \\
# \end{pmatrix}
# \end{align}

# <codecell>

s0=''.join(['G', 'A', 'G', 'A', 'C', 'T', 'C', 'A', 'T', 'A', 'A', 'T', 'G',
        'T', 'C'])
s1=''.join(['G', 'T', 'T', 'A', 'G', 'A', 'C', 'A', 'A', 'T', 'A', 'T', 'G',
        'T', 'T'])
s2=''.join(['C', 'T', 'G', 'A', 'A', 'A', 'C', 'A', 'C', 'T', 'A', 'T', 'G',
        'T', 'C'])
s3=''.join(['C', 'T', 'G', 'A', 'C', 'A', 'C', 'G', 'T', 'T', 'A', 'T', 'G',
        'T', 'C'])
s4=''.join(['G', 'G', 'G', 'A', 'C', 'A', 'C', 'A', 'T', 'A', 'A', 'T', 'G',
        'T', 'C'])
s5=''.join(['G', 'T', 'G', 'A', 'G', 'A', 'G', 'A', 'C', 'G', 'A', 'T', 'G',
        'T', 'C'])
s6=''.join(['G', 'T', 'G', 'A', 'C', 'A', 'C', 'A', 'T', 'T', 'A', 'T', 'G',
        'T', 'C'])
s7=''.join(['G', 'T', 'G', 'A', 'T', 'A', 'T', 'G', 'T', 'T', 'A', 'T', 'G',
        'A', 'C'])
s8=''.join(['G', 'T', 'G', 'A', 'C', 'A', 'C', 'A', 'T', 'G', 'G', 'G', 'G',
        'T', 'C'])
s9=''.join(['G', 'T', 'G', 'A', 'C', 'A', 'C', 'A', 'T', 'T', 'A', 'T', 'G',
        'T', 'C'])
print s0,s1,s2,s3,s4
print s5,s6,s7,s8,s9

# <markdowncell>

# Vypísanie nájdených sekvencií z fasta súboru. Treba na otestovanie, či sa naozaj nachádzajú v súbore

# <headingcell level=2>

# Hlavné vykonanie gibbsovho vzorkovania pre MotSamShort

# <markdowncell>

# Pustanie gibbsovho vzorkovania na MotSamShort.fasta s poctom baz 9 a maximalne 4 mutaciami. Do premennych motif,score,array ukladame navratove hodnoty z metody.Tie zodpovedaju consensus motivu, skore consensus motivu od ostatnych sekvencii => sucet vzdialenosti od ostatnych sekvencii a vypis najdenych sekvencii. 
# 
# Kedze skore obsahuje sucet vzdialenosti consensus motivu od ostatnych sekvencii a vieme ze sa maju vyskytovat najviac 4 mutacie, tak jedini vhodni kandidati su taki, kde skore bude menej ako 40 (pocetsekvencii * pocetmutacii). Kedze ale porovnavame menej baz ako v predchodzom pripade, tak dost casto dostaneme riesenie, v ktorom bude skore mensie ako 40, dokonca mensie ako 30. Tych najlepsich kandidatov si ukladame do 10 prvkoveho pola aby sme mohli vidiet ako sa vyvijalo hladanie riesenia.

# <codecell>

# Nacitanie sekvencii z fasta suboru
sequences2 = readFasta(seq2file,verbose=0);
# ziskanie hodnot
values = sequences2.values()
# a klucov, ktore nevyuzijeme
keys = sequences2.keys()

# Vytvorenie pola do ktoreho budeme ukladat 10 najlepsich behov
best2=[['',150,[]] for i in range(10)]

timeArray=[]

# 50 krat opakujeme gibbsovo samplovanie s 
for j in range(50) :
    # Nastavenie casu
    start = time.clock()
    # Samplovanie pustame s dlzkou motivu 15. Argument 200 znaci kolko krat opakujeme iterovanie motivmi
    motif,score,array = gibbsSampling(values,9,4,200,choice=0)
    # Vypisanie aktualnej iteracie opakovania a ziskaneho skore
    print 'Iter:',j,'with score',score,'and motify',motif
    # Nasledne priradim tuto trojicu do pola best na index kam patri. Ked pole pretecie tak vymazem posledny prvok (ten najvacsi)
    for k in range(10) :
        if score < best2[k][1] :
            best2.insert(k,[motif,score,array])
            del best2[-1]
            break
    timeArray.append(time.clock()-start)        

print 'Gibbsovo vzorkovanie zbehlo v case',sum(timeArray)             

# <headingcell level=3>

# Vysledky pre MatSamShort

# <codecell>

# Tento skript vyberie kazdu trojicu v ktorej su ulozene najlepsie najdene motify.
# Z kazdej trojice vyberie pole (na druhom mieste) v ktorych su kusky sekvencii a 
# pre kazde zistime hamming vzdialenost od consensus motifu. Ked su vzdialenosti od consensus motifu
# mensie alebo rovne hodnote mut tak je to spravny motif

mut = 4

for i in range(len(best2)) :
    w=[hamDist(best2[i][0],seq) for seq in best2[i][2]]    
    if all([x <= mut for x in w]) :
        print 'Optimal',w,'Consensus motif',best2[i][0]
    else :
        print w,'Consensus motif',best2[i][0]

# <markdowncell>

# Vypocet bezi rychlejsie ako v predchadzajucom pripade hlavne kvoli zvoleniu poctu bazi na 9.
# 
# Vyššie je vidieť 10 najlepších kandidátov. V tomto pripade tiez nemozme vzdy na 100% povedat ci prvý (dokonca ani ostatok) z nich je ten optimálny. V tomto pripade prve riesenie (ktore ma najmensi sucet vzdialenosti od consensu motivu) ma vzdialenosť consensus motivu od najdenych sekvencii [4, 2, 1, 2, 1, 2, 2, 3, 1, 2], kde su najviac 4 mutacie. Druhy, treti, stvrty, siesty, deviaty pripad je to iste, no ako vidiet tak dostavame vzdy rozne consensus motivy.
# 
# Piaty, siedmy, osmy a desiaty pripad nie su v tomto pripade tie optimalne aj ked napriklad 5 riesenie ma mensi sucet vzdialenosti ako nasledujuce.

# <codecell>

best2

# <headingcell level=2>

# Konecne Vysledky

# <markdowncell>

# 1.
# V subore MotSam sme nasli motiv GTGACACATTATGTC. 
# 
# Parametre gibbsovho vzorkovania boli : 
# 
#     pocet opakovani = 50
#     
#     pocet iteracii = 200
#     
#     pocet bazi = 15
#     
#     pocet mutacii = 4
# 
# 2.
# Vypocet v tomto pripade trval priblizne 3 minuty. V priebehu vypoctu bolo po kazdom zbehnuti gibbsovho vzorkovanie vidiet ten najlepsi najdeny motiv s vypocitanym skore pre najdene sekvencie. Skore sa pohybovali do 60 a najlepsi motiv [GTGACACATTATGTC] sa zda byt so skore 25. Neukazujem startovacie pozicie sekvencii, kedze to je len formalita.
# 
# 3.
# Nasli sme teda motiv GTGACACATTATGTC, so skore 25
# 
# Vzdialenost tohoto motivu od ostatnych sekvencii je [3, 4, 3, 2, 2, 4, 0, 4, 3, 0]. (sucet dava 25)
# 
# 4.
# Gibbsovo vzorkovanie ma mnoho problemov. 
# 
# Problemy :
# 
# Vo vela pripadoch bezi dobre, ale velmi lahko dokaze skonvergovat k nejakemu suboptimalnemu rieseniu. (nespolahlivost)
# 
# Dalsim problemom je, ze algoritmus je pritahovany k low-complexity oblastiam hlavne v pripade, ked vstupne sekvencie obsahuju niektore nucleotidy viacero razy nez ine.
# 
# Riesenia :
# 
# Opakovanie gibbsovho vzorkovania viacero razy a dufat ze niektory z behov da optimalne riesenie. Toto sme robili aj my (pocetopakovanie=50, spolahlivejsie). Taktiez vyuzivat Laplaceovo pravidlo uspechu, kde pri vytvarani profilu (matice) pridame ku kazdej bunke jednotku => nedostaneme pravdepodobnost 0 v ziadnom pripade.
# 
# Vyuzivat relativne entropie namiesto frekvencii (zabranenie druheho pripadu)
# 
# Rychlost by sa dala vyriesit znizenim poctu iteracii v jednom gibbsovom vzorkovani, ci znizenim poctu opakovani vzorkovania. Tak ci tak znizime kvalitu najdeneho riesenia.
# 
# 5.
# V subore MotSamShort sme nasli viacero motivov, ktore pasuju podmienke aby boli najviac 4 mutacie.
# 
# motivy = [TTAGAGATT,CTTCGCGCA,TAGAGCATT,CTTGGTATT,GTGGTTCCA,AGCGTACGT]
# 
# Najmensie skore ma znich motiv TTAGAGATT, no to neznamena, ze to ma byt ten implantovany.
# 
# Vypocet v tomto pripade bol rychlejsi (90 sekund), pri ponechani rovnakych parametrov ako v prvom pripade. Vysledok je ale nejednoznacny kvoli viacero najdenym vyhovujucim motivom. Taktiez sa mohlo stat ze skutocny motiv dokonca neni zaradeny medzi kandidatmi (10 prvkove pole najlepsich kandidatov), kvoli tomu ze bol vytlaceny nejakym s mensim skore.
# 
#     

