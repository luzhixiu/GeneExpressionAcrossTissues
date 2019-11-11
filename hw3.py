#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 02:58:12 2019

@author: lu
"""

import numpy as np
import random
import math


import matplotlib.pyplot as plt
t=np.linspace(-10,10,100) # 100 linearly spaced numbers
y=[]
for i in t:
    j= 50* math.sin(i)+ (i*i)
    y.append(j)
plt.plot(t,y)

#c = 0.1, 0.2, 1; epsilon = 0.01, 0.1
def gd(x,c,epsilon):
    # c is learning rate, epsilon controls accuracy
    result=[]
    while(True):
        result.append(x)
        x1= x-(((50*math.cos(x)+2*x))*c)
        if abs(x1-x)<epsilon:
            break
        else:
            x=x1
    return result

x_inrange=[]
result=gd(1,0.1,0.01)
for i in result:
    if abs(i)<=10:
        x_inrange.append(i)
#        
y=[]
for i in x_inrange:
    j= 50* math.sin(i)+ (i*i)
    y.append(j)

for i in range(len(y)):
    print x_inrange[i],
    print y[i]
    
print np.min(y)



