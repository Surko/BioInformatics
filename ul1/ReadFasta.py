# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

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

D=readFasta('d:\\cvicenia\\BioInformatics\\ul1\\test.fasta',0)

# <codecell>


