#!/usr/bin/python2.5
#-*- coding:utf-8 -*-

import math, random, numpy as np


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
    Vp = np.zeros(V.shape)
    for idx in range(len(V)):
        Vp[idx] = discreto(V[idx],mx,mn,n)
    return (Vp)


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
    N = len(X)
    M = []

    #contamos las ocurrencias
    for x in X:
        D[x] = D.get(x,0)+1

    for d in D:
        D[d] = D[d] / N

    #determinamos probabilidades
    M = list(filter( lambda x: D[x]>0, D))

    return D


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
    PX = prob(X)
    PY = prob(Y)

    D = dict()
    tot = 0

    #aqui contamos los simbolos
    for x,y in zip(X,Y):
        D[(x,y)] = D.get( (x,y), 0 ) + 1
        tot += 1

    #determinamos la probabilidad de aparicion de cada
    #elemento conjunto
    for idx in D:
        D[idx] /= float(tot)

    return [PX, PY, D] 

def H(X,base=2):
    """
    Regresa la Entropia de Shannon de p en los estados que deben venir
    en un diccionario como esta definido in prob, las unidades de
    H esta en bits

    regresa H: float con la entropia de X
    """
    H = sum( -1* np.log(np.array(list(X.values())))/np.log(base))

    return H

def Hxy(P,base=2):
    """
    calcula la entropia conjunta
    recibe P: entropia de estados conjuntos calculada en prconj

    regresa la entropía conjunta entre dos sistemas

    """
    return H(P,base)

def entropia(X):
    """
    Calcula la entropia de la secuencia en el atributo X
    esta puede ser una cadena
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
    PX, PY, PXY = prconj(X,Y)
    #print(PX)
    h1, h2 = H(PX), H(PY)
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

def mapea_val( S, el ):
    """
    Asigna la entropia puntual de el respecto
    de los valores en S donde se encuentra
    contenido
    """
    pass

