#!-*-coding:utf8-*-
import sys
sys.path.append('~/checo/lib/')
import numpy as np
import inf_mutua as im

T = np.arange(0,4*np.pi,0.01)
Y = np.sin(T)
Yd = im.discretiza(Y,10)
Yv = Yd[0:100]
I = []
for i in range(len(Yd)-100):
    V = Yd[i:i+100]
    im = im.im(V,Yv)
    I.append(im)


