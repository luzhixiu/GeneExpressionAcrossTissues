#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 04:58:51 2019

@author: lu
"""

#% learning process (training process)
#finish = 0;
#while ~finish
#disp(w)
#for i=1:tr
#z(i) = w(1:input) * x(i,:)' > w(input+1);
#w(1:input) = w(1:input) + (T(i)-z(i))* x(i,:);
#w(input+1) = w(input+1) - (T(i)-z(i));
#end
#if sum(abs(z-T')) == 0
#finish = 1;
#end
#disp(z)
#pause
#end
#disp(w)
import numpy as np
import random as rd
test="off" #switch this to on to see out put of iterations
input_=2
tr=4
w=rd.randint(1,input_+1)
x=np.array([[0,0],[1,0],[0,1],[1,1]])
AND=np.array([0,0,0,1]) #expected output for AND
OR=np.array([0,1,1,1])  #expected output for OR
XOR=np.array([0,1,1,0]) #expected output for XOR

import numpy as np

class Perceptron(object):

    def __init__(self, no_of_inputs, iterations=10):
        self.iterations = iterations
        self.weights = np.zeros(no_of_inputs + 1)
           
    def predict(self, inputs):
        sum = np.dot(inputs, self.weights[1:]) + self.weights[0]
        if sum > 0:
          activation = 1
        else:
          activation = 0  
        global test
        if "on" in test:
            print activation          
        return activation

    def train(self, training_inputs, labels):
        for _ in range(self.iterations):
            for inputs, label in zip(training_inputs, labels):
                prediction = self.predict(inputs)
                self.weights[1:] += (label - prediction) * inputs
                self.weights[0] += (label - prediction)
                

perceptron = Perceptron(input_)

print "Test inputs are: [0,0] [0,1] [1,0],[1,1]"
print "Result for AND"
perceptron.train(x, AND)
#print perceptron.weights
print perceptron.predict([0,0])
print perceptron.predict([0,1])
print perceptron.predict([1,0])
print perceptron.predict([1,1])


print "Result for OR"
perceptron.train(x, OR)
print perceptron.predict([0,0])
print perceptron.predict([0,1])
print perceptron.predict([1,0])
print perceptron.predict([1,1])

print "Result for XOR"
perceptron.train(x, XOR)
print perceptron.predict([0,0])
print perceptron.predict([0,1])
print perceptron.predict([1,0])
print perceptron.predict([1,1])
