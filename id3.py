#!/usr/bin/python
#!-*-coding:utf8-*-

from numpy import log2

def H(X):
    """
    Regresa la entropia de un conjunto de frecuencias
    relativas 
    """
    S = 0
    for x in X:
        if x>0:
            S += -x*log2(x)
        else:
            S += 0
    return S

def entropia(V):
    """
    Regresa la entropia del conjunto V
    donde primero cuenta las frecuencias
    de cada valor v
    """
    S = list(V)
    X = dict()
    for v in S:
        X[v] = X.get(v,0)+1
    L = [ X[x] for x in X ]
    L = [float(x)/sum(L) for x in L]
    return H(L)

def ganancia(C,iA,iS):
    """
    Determina la ganancia sobre S=C[:,iS] del atributo A=C[:,iA]
    """
    S = C[:,iS]
    tmS = float(len(S))
    A = C[:,iA]
    G=0
    V= set(A)
    for v in V:
        idxs = C[:,iA]==v
        sbc = C[idxs]
        tm = len(sbc)
        G += tm/tmS * entropia(sbc[:,iS])
    return entropia(S) - G
