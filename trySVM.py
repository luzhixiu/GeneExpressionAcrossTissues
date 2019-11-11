#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 06:14:24 2019

@author: lu
"""


from sklearn.svm import SVC
import numpy as np
X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
Y = np.array([1, 1, 2, 3])
clf = SVC(gamma='auto')
clf.fit(X,Y)
print clf.predict([[2,1]])