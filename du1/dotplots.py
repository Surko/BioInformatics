# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

vectorx='MAAPSRTTLMPPPFRLQLRLLILPILLLLRHDAVHAEPYSG'
vectory=vectorx
line=[];
for i in range(len(vectory)) :
    for j in range(len(vectorx)) :
        if vectory[i]==vectorx[j] :
            line.append('■');
        else :
            line.append(' ');
    print ''.join(line)
    line =[];

# <codecell>

from numpy import *
import Image as img
import ImageOps
import ImageFilter
import random

def smscale(input, mn, mx,inverse = 0) :    
    if inverse==1 :
        input = [-x for x in input]
    maxElem = max(input)
    minElem = min(input)
    q = (mx - mn)/(maxElem - minElem)
    for i in range(len(input)) :
        input[i] = mn + (input[i]-minElem)*q
    return input

def windowFilled(win1,win2) :
    common = 0;
    for i in range(len(win1)) :
        for j in range(len(win2)) :
            if win1[i]==win2[j] :
                common = common + 1
    return common

def windowFilledDiagonal(win1,win2) :
    common = 0;
    for i in range(min(len(win1),len(win2))) :
        if win1[i]==win2[i] :
            common = common + 1
    return common

def dotPlot(inputVector,minGray=0.0,maxGray=1.0,background=0.0,n=1,threshold=0,diagonal=0,distribution='uniform',scaling=0,choice=0) :
    inputLength = len(inputVector);
    
    if distribution =='uniform' :
        if diagonal == 0 :
            distribution = [i * (maxGray - minGray) / n / n for i in range(n * n + 1)]
        else :
            distribution = [i * (maxGray - minGray) / n  for i in range(n + 1)]
    elif distribution == 'beta' :        
        if diagonal == 0 :
            distribution = [ (maxGray - minGray) * random.betavariate(1,0.5) for i in range(n * n + 1)]
        else :    
            distribution = [ (maxGray - minGray) * random.betavariate(1,0.5) for i in range(n + 1)]
        distribution.sort()
        if (scaling == 1) :
            distribution = smscale(distribution,minGray,maxGray)
    else :
        distribution = [ (maxGray - minGray) * x for x in distribution ]        
        distribution.sort()
        if (scaling == 1) :
            distribution = smscale(distribution,minGray,maxGray)        
    
    raster = {}    
    if diagonal==0 :        
        for i in range(n*n+1) :
            raster[i]=(int)((minGray+distribution[i])* 255)    
    else :        
        for i in range(n+1) :
            raster[i]=(int)((minGray+distribution[i])* 255)        
    
    print raster
    
    arr = ones((inputLength,inputLength)) * int(background * 255)
    
    for i in range(inputLength) :
        if i - n/2 < 0 :
            win1=inputVector[:i+1+n/2]
        elif i + n/2 > inputLength :
            win1 = inputVector[i-n/2:]
        else :
            win1 = inputVector[i-n/2:i+1+n/2]
        for j in range(inputLength) :
            if j - n/2 < 0 :
                win2=inputVector[:j+n/2]
            elif j + n/2 > inputLength :
                win2 = inputVector[j-n/2:]
            else :
                win2 = inputVector[j-n/2:j+1+n/2]
            if diagonal==0 :
                count = windowFilled(win1,win2)
            else : 
                count = windowFilledDiagonal(win1,win2)
            if count > threshold :
                if (choice == 1) :
                    arr[i,j]=raster[(count - threshold) * n/(n - threshold)] 
                else :
                    arr[i,j]=raster[count - threshold] 
    im = img.fromarray(arr)
    im.show()
    
vectory = 'GCTGAGCTTACGAACAAAACTGCAGCAGTGATGTAAATATACGTACGGTATTTTAGGCTTGTACACCCCTCTATTACACATACACACACGCACACACACACACACACACCTGAGGTTACTGAAGTAAGGTTGGAGACGGTACTTGTCTATCTCCCAGCCGAAGTGGTCTTCCGCTGAGCAGAGTTCCTTTGCCACCCTGAATCATGGCTGTTGGTTCAGTATGAAGGTGTTGTACCCAGTAGCGCTGTTTTACAGACACACACACAAACGGGCACACACACATATACACACAGATACACACACACACACTAGGTAACTTATTCTTTGAACTTTTCTATCT'
vectorx = 'MAAPSRTTLMPPPFRLQLRLLILPILLLLRHDAVHAEPYSGGFGSSAVSSGGLGSVGIHIPGGGVGVITEARCPRVCSCT GLNVDCSHRGLTSVPRKISADVERLELQGNNLTVIYETDFQRLTKLRMLQLTDNQIHTIERNSFQDLVSLERLDISNNVI TTVGRRVFKGAQSLRSLQLDNNQITCLDEHAFKGLVELEILTLNNNNLTSLPHNIFGGLGRLRALRLSDNPFACDCHLSW LSRFLRSATRLAPYTRCQSPSQLKGQNVADLHDQEFKCSGLTEHAPMECGAENSCPHPCRCADGIVDCREKSLTSVPVTL PDDTTDVRLEQNFITELPPKSFSSFRRLRRIDLSNNNISRIAHDALSGLKQLTTLVLYGNKIKDLPSGVFKGLGSLRLLL LNANEISCIRKDAFRDLHSLSLLSLYDNNIQSLANGTFDAMKSMKTVHLAKNPFICDCNLRWLADYLHKNPIETSGARCE SPKRMHRRRIESLREEKFKCSWGELRMKLSGECRMDSDCPAMCHCEGTTVDCTGRRLKEIPRDIPLHTTELLLNDNELGR ISSDGLFGRLPHLVKLELKRNQLTGIEPNAFEGASHIQELQLGENKIKEISNKMFLGLHQLKTLNLYDNQISCVMPGSFE HLNSLTSLNLASNPFNCNCHLAWFAECVRKKSLNGGAARCGAPSKVRDVQIKDLPHSEFKCSSENSEGCLGDGYCPPSCT CTGTVVACSRNQLKEIPRGIPAETSELYLESNEIEQIHYERIRHLRSLTRLDLSNNQITILSNYTFANLTKLSTLIISYN KLQCLQRHALSGLNNLRVVSLHGNRISMLPEGSFEDLKSLTHIALGSNPLYCDCGLKWFSDWIKLDYVEPGIARCAEPEQ MKDKLILSTPSSSFVCRGRVRNDILAKCNACFEQPCQNQAQCVALPQREYQCLCQPGYHGKHCEFMIDACYGNPCRNNAT CTVLEEGRFSCQCAPGYTGARCETNIDDCLGEIKCQNNATCIDGVESYKCECQPGFSGEFCDTKIQFCSPEFNPCANGAK CMDHFTHYSCDCQAGFHGTNCTDNIDDCQNHMCQNGGTCVDGINDYQCRCPDDYTGKYCEGHNMISMMYPQTSPCQNHEC KHGVCFQPNAQGSDYLCRCHPGYTGKWCEYLTSISFVHNNSFVELEPLRTRPEANVTIVFSSAEQNGILMYDGQDAHLAV ELFNGRIRVSYDVGNHPVSTMYSFEMVADGKYHAVELLAIKKNFTLRVDRGLARSIINEGSNDYLKLTTPMFLGGLPVDP AQQAYKNWQIRNLTSFKGCMKEVWINHKLVDFGNAQRQQKITPGCALLEGEQQEEEDDEQDFMDETPHIKEEPVDPCLEN KCRRGSRCVPNSNARDGYQCKCKHGQRGRYCDQGEGSTEPPTVTAASTCRKEQVREYYTENDCRSRQPLKYAKCVGGCGN QCCAAKIVRRRKVRMVCSNNRKYIKNLDIVRKCGCTKKCY';
dotPlot(vectory,0.3,1.0,n=17,threshold=7,diagonal=1,distribution='uniform',choice=0)

# <codecell>

import Image as img
vectorx = 'MAAPSRTTLMPPPFRLQLRLLILPILLLLRHDAVHAEPYSGGFGSSAVSSGGLGSVGIHIPGGGVGVITEARCPRVCSCT';

mat = [255 *int(x==y) for x in vectorx for y in vectorx]
im = img.new("L",(len(vectorx),len(vectorx)))
im.putdata(mat)
im.show()

# <codecell>


