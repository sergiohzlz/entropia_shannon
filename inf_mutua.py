#!/usr/bin/python2.5
#-*- coding:utf-8 -*-

import math, random


def discreto(var, vmax, vmin, n_edos):
    """
    discretiza la variable v dentro
    del rango establecido por max y  min 
    con la cantidad de estados definidos en n_edos

    regresa dig: un entero que var discretizada
    """
    intervalo = vmax - vmin
    R = float(intervalo) / float(n_edos)

    dig = int((var-vmin) / R) - 1
    if dig < 0:
        return 0

    return dig


def discretiza(V, n): 
    """
    Esta función los valores descritos en V en 
    n estados distintos
    """
    mx = max(V)
    mn = min(V)
    rango = mx - mn
    dt = 1.*rango / n
    Vp = V + mn
    Vp /= dt
    return int(Vp)


def prob(X):
    """
    determina la probabiliad de los elementos en 
    el conjunto X

    regresa [P, N, M]: un vector con 3 objetos:
    P. probabilidad de ocurrencia de cada simbolo en X
    N. la cantidad original de elementos en X
    M. cantidad de simbolos con probabilidad no nula
    """
    D = {} #diccionario donde guardamos la ocurrencia de cada simbolo
    P = {} #diccionario donde van a estar las probabilidades
    N = len(X)
    M = []

    #contamos las ocurrencias
    for x in X:
        D[x] = D.get(x,0)+1

    #determinamos probabilidades
    for i in range(len(D)):
        P[D.keys()[i]] = 0
        P[D.keys()[i]] = float(D.values()[i]) / N

    M = filter( lambda x: P[x]>0, P)

    return [P, N, len(M)]


def prconj(X,Y):
    """
    determina la probabilidad conjunta 
    entre los conjuntos X e Y que deben estar
    discretizados

    Calcula 
    p(X,Y)=p(X=x,Y=y)

    regresa [PX, PY, PXY]:
    PX: probabilidad de aparicion de X
    PY: probabilidad de aparicion de Y
    PXY: probabilidad conjunta de los sistemas X e Y
    """
    [PX,NX,MX] = prob(X)
    [PY,NY,MY] = prob(Y)

    D = dict()
    tot = 0

    #aqui contamos los simbolos
    for i in range(NX):
        for j in range(NY):
            D[X[i],Y[i]] =  D.get((X[i],Y[i]),0) + 1
            tot += 1

    #determinamos la probabilidad de aparicion de cada 
    #elemento conjunto
    for idx in D:
        D[idx] /= float(tot)

    #prob no nula
    M = filter( lambda x: D[x]>0, D)

    return [PX, PY, D, len(M)]

def H(X):
    """
    Regresa la Entropia de Shannon de p en los estados que deben venir 
    en un diccionario como esta definido in prob, las unidades de
    H esta en bits

    regresa H: float con la entropia de X
    """
    H = 0.0

    for j in range(len(X)):
        logP = 0.0;
        if (X.values()[j] > 0):
            logP = math.log(X.values()[j],2)
        H += X.values()[j] * logP
    #H = -H / log(2)

    return -H

def Hxy(P):
    """
    calcula la entropia conjunta
    recibe P: entropia de estados conjuntos calculada en prconj

    regresa la entropía conjunta entre dos sistemas
    """
    HXY = 0.0
    for idx in P:
            logP = 0.0;
            if P[idx] > 0:
                    logP = math.log(P[idx],2)
            HXY += P[idx] * logP
    #HXY = -HXY / math.log(2)

    return -HXY


def Hmarginal(X):
    """
    Calcula la entropia marginal del conjunto de datos X
    """
    P = prob(X)
    return H(P[0])

def Hconj(X,Y):
    """
    calcula la entropia conjunta de los conjuntos de datos
    X e Y
    """
    P=prconj(X,Y)
    return Hxy(P[2])


def inf_mutua(H1, H2, HXY):
    """
    calcula la entropia conjunta de H1 y H2
    recibe tambien HXY que es la entropia conjunta de H1 y H2

    regresa la informacion mutua dependiendo de la entropia
    del primer sistema, la del segundo y la conjunta
    """
    I = H1 + H2 - HXY
    return I

def im(X,Y):
    """
    calcula la información mutua entre el conjunto de datos 
    X e Y 
    """
    h1 = Hmarginal(X)
    h2 = Hmarginal(Y)
    hxy = Hconj(X,Y)
    return inf_mutua(h1,h2,hxy)

def imc(X,Y):
    """
    Calcula la correcion de informacion mutua con el metodo
    de elemento finito
    """
    PX = prob(X)
    PY = prob(Y)
    PXY = prconj(X,Y)
    N = PX[1]
    delta = (float(PXY[3]-PX[2]-PY[2]+1))/(2*N)
    I = im(X,Y)
    Ic = I + delta
    return (Ic, I, delta)

def ganancia(S,A,verbose=False):
    """
    calcula la funcion Ganancia con el conjunto de datos de la clase S
    y el conjunto de atributos A
    """
    hs = Hmarginal(S)
    P=prob(A)
    G=dict()
    for val in zip(A,S):
        L=G.get(val[0],[])
        L.append(val[1])
        G[val[0]]=L
    if verbose: print(P[0]); print(G)
    H=0.0
    for val in P[0]:
        if verbose: print(val); print(P[0][val]); print(G[val]); print(Hmarginal(G[val])); print(P[0][val]*Hmarginal(G[val]))
        H += P[0][val]*Hmarginal(G[val])

    return hs-H

