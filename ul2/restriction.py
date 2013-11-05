# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from itertools import combinations as comb
def formDelta(array) :
    length = len(array)
    deltaX=[array[j+i+1]-array[i] for i in range(length) for j in range(length-i-1)]
    deltaX.sort()
    return deltaX

def anotherBFPDP(L,n) :
    maximum = max(L)    
    fullL=[0]+L
    
    for xList in comb(fullL,n) :
        if L==formDelta(xList) :
            return xList
    return []

def subset(X,Y) :
    return all([X.count(c)<=Y.count(c) for c in set(X)])

def dFunction(y,X) :
    return [abs(y-x) for x in X]

def place(L,X,W):
    if L==[] :
        X.sort()
        print X
        return
    maximum=max(L)
    dElem = dFunction(maximum,X)
    if subset(dElem,L) :
        X=X+[maximum]
        temp=L
        for c in dElem :
            L.remove(c)
        place(L,X,W)
        X.remove(maximum)        
        L=temp
    dElem = dFunction(W-maximum,X)
    if subset(dElem,L) :
        X=X+[W-maximum]
        temp=L
        for c in dElem :
            L.remove(c)
        place(L,X,W)
        X.remove(W-maximum)        
        L=temp
    return

def partialDigest(L) :    
    width = max(L)
    L.remove(width)
    X=[0,width] 
    place(L,X,width)    
        
def SPDP(full,shorter):
    full.sort()        
    if (len(full) % 2) != 0 : 
        error(' Zle velkosti ')
    n = len(full)/2
    X=(n)*[0]    
    suma=full[0]+full[2*n-1]
    compShorter = (n+1) * [0]
    order=2*n*[0]    
    for i in range(n) :
        order[i]=(full[i],full[2*n-i-1])
        order[2*n-i-1]=(full[2*n-i-1],full[i])
        
    for k in comb(range(2*n),n) :        
        combination=n*[0]
        for i in range(n) :
            X[i]=order[k[i]][0]                        
            combination[i]=order[k[i]]
                                 
        compShorter[0]=X[0]
        for i in range(n-1) :
            compShorter[i+1]=X[i+1]-X[i]
        compShorter[n]=suma-X[n-1]        
        if compShorter == shorter :
            print X,'with orderings',combination,'with indexes',k
    
    
    

# <codecell>

anotherBFPDP([2,998,1000],3)

# <codecell>

partialDigest([2, 2, 3, 3, 4, 5, 6, 7, 8, 10])

# <codecell>

SPDP([2,14,3,13,7,8,8,9],[2,6,1,4,3])

# <codecell>


# <codecell>


