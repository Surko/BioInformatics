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

# <codecell>

anotherBFPDP([2,998,1000],3)

# <codecell>

partialDigest([2, 2, 3, 3, 4, 5, 6, 7, 8, 10])

# <codecell>

a=[]
a==[]

# <codecell>


