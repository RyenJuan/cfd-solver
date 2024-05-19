#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tues Apr 19 17:39:00 2022

@author: ryanhuang
"""

import matplotlib.pyplot as plt 
import numpy as np
from numpy.linalg import inv


# number of gridpoints
N = 100
    
# boundary conditions
T1, T2 = 0, 0

# time step and space step
DT, DX = 0.00001, 0.01

r = DT / (DX**2)   # r = 0.1

# constants
T_i     = 1
B_i     = 10


def heat_eq(T1, T2, T3, non_homo):
    '''
    Discretized heat equation to calculate the next temperature at a specific grid point
    '''
    if non_homo:
        new_T = r * T1 + (1 - 2*r) * T2 + r * T3
    else:
        new_T = r * T1 + (1 - 2*r) * T2 + r * T3 - 2 * DT
    
    return new_T


def update(temperatures, non_homo):
    ''' 
    Input a list of temperature values
    
    Output a new list of temperature values for the next timestep based on the discretized heat equation
    
    '''
    temp = []
    
    for i in range(1, len(temperatures) - 1):
        if i - 1 == 0: 
            temp.append(heat_eq(T1, temperatures[i], temperatures[i+1], non_homo))
        elif i + 1 == 100:
            temp.append(heat_eq(temperatures[0], temperatures[i], T2, non_homo))
        else:
            temp.append(heat_eq(temperatures[i-1], temperatures[i], temperatures[i+1], non_homo))
    
    return temp


def main():
    '''
    --------------------------------------------------------------------------
    implicit method
    
    '''
    

        
    # setup the RHS and LHS matrices
    a = np.zeros((N+1, N+1))
    b = np.zeros((N+1,1))
    b.fill(T_i)
    
    for i in range(0, N-1):
        a[i+1][i] = -r
        a[i+1][i+1] = (1 + 2*r)
        a[i+1][i+2] = -r
        
    # boundary conditions
    a[0][0] = 1
    a[0][1] = -1
    a[1][0] = 0
    a[1][1] = 1 + r
    
    a[N][N-1] = 1
    a[N][N] = 1 + DX*B_i 

    a_inv = inv(a)

    for i in range(5000):
        new = np.matmul(a_inv, b)
        b = new
        b[0] = b[1]

    
    plt.figure(1)
    plt.title(f"Implicit Method at B_i = {B_i}")
    plt.xlabel("Points")
    plt.ylabel("Dimensionless Temperature")
    plt.plot(list(range(len(b.transpose()[0]))), b.transpose()[0])
    plt.axis([0, N, 0, 1.5])
    plt.show()

    
    
    '''
    --------------------------------------------------------------------------
    explicit method
    
    '''
    
    # initialize temperature values at T_i
    T = [T1] + [1 for i in range(N - 2)] + [T2]

    # run for 5000 steps
    for i in range(5000):
        
        new = [T1] + update(T, non_homo = True) + [T2]
        T = new
        T[0] = T[1]
        
        T[N-1] = T[N-2]/(1+0.01*DX)
        
    plt.figure(2)
    plt.title(f"Explicit Method at B_i = {B_i}")
    
    # plot the final temperature distribution
    plt.plot(list(range(1,N+1)), list(reversed(T)))
    plt.axis([0, N, 0, 1.5])
    plt.xlabel("Points")
    plt.ylabel("Dimensionless Temperature")
    plt.show()

        

if __name__ == "__main__":
    main()











