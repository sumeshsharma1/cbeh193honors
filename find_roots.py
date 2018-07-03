import numpy as NP
import itertools

#Instructions

#Example for inputting the names:
# input ['hermite','legendre']
# hermite for gaussian variables, legendre for uniform variables 

#order = int(input("Order? "))
order = 3
#np = int(input("Number of parameters "))
np = 2
#names = input("Probability descriptions array ")
names = ['hermite','hermite']

#define coefficients for polynomial of degree 1 higher than order
coeff = NP.zeros(order+1); coeff = NP.append(coeff, 1.0)
# generate roots of for a polynoimal of degree 1 higher than order
roots = []
for i in range(np):
    if names[i] == 'hermite':
        roots += [NP.polynomial.hermite_e.hermeroots(coeff).tolist()]
    elif names[i] == 'legendre':
        roots += [NP.polynomial.legendre.legroots(coeff).tolist()]
# construct sample point set
ns = pow((order+1),np) # number of samples
sample_points = NP.zeros((ns,np)) # set of sample points
j = 0
for point in itertools.product(*roots):
    sample_points[j,:] = NP.asarray(point)
    j += 1

print(sample_points)

