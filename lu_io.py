#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 10:21:47 2019

@author: lu
"""

def writeListToCsv(ls,fileName):
    f=open(fileName,"w+")
    for item in ls:
        f.write(str(item))
        f.write(",")
        f.write("\n")
        
    f.close()
    