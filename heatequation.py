#!/usr/bin/env python3
"""
Created on Thu Jan 20 14:06:57 2022

@author: ryanhuang
"""
import matplotlib.pyplot as plt 
import numpy as np
from numpy.linalg import inv


#number of gridpoints
N = 100
    
# boundary conditions
T1, T2 = 0, 0

# time step and space step
DT, DX = 0.0000001, 0.001

r = DT / (DX**2)   # r = 0.1

# r=0.501


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
        # print(i)
        if i - 1 == 0: 
            temp.append(heat_eq(T1, temperatures[i], temperatures[i+1], non_homo))
        elif i + 1 == 100:
            temp.append(heat_eq(temperatures[0], temperatures[i], T2, non_homo))
        else:
            temp.append(heat_eq(temperatures[i-1], temperatures[i], temperatures[i+1], non_homo))
    
    return temp


def main():
    '''
    -----------------------------------------------------------------------------------------
    part a: steady state
    
    '''
##    
##    a = np.zeros((N+1, N+1))
##    b = np.zeros((N+1,1))
##    
##    a[0][0] = 1
##    a[N][N] = 1
##    
##    b[N][0] = 1
##    
##    
##    for i in range(0, N-1):
##        a[i+1][i] = 1
##        a[i+1][i+1] = -2
##        a[i+1][i+2] = 1
##        
##    a_inv = inv(a)
##    
##    result = np.matmul(a_inv, b)
##    
##    result_transpose = np.transpose(result)
##    
##    plt.figure(1)
##    plt.title("Steady State Approach")
##
##    plt.plot(list(range(len(result_transpose[0]))), result_transpose[0])
##    plt.plot(result_transpose)
##            
##    
##    
##    '''
##    -----------------------------------------------------------------------------------------
##    part b: transient
##    
##    '''
##    
##    # initialize temperature values at 0
##    T = [T1] + [1 for i in range(N - 2)] + [T2]
##
##    # run for 100000 steps
##    for i in range(5000):
##        
##        new = [T1] + update(T, non_homo = True) + [T2]
##        T = new
##        T[0] = T[1]
##        
##        T[N-1] = T[N-2]/(1+0.01*DX)
##        
##        plt.figure(2)
##        plt.title("Explicit Method")
##        
##        # plot the final temperature distribution
##        plt.plot(list(range(1,N+1)), list(reversed(T)))
##        plt.axis([0, N, 0, 1.5])
##        plt.show()
##        plt.pause(0.01)
    
    
        
    '''
    -----------------------------------------------------------------------------------------
    part c: non-homogeneous steady state
    
    '''
    
    a = np.zeros((N+1, N+1))
    b = np.zeros((N+1,1))
    
    b.fill(2*(DX**2))
    
    a[0][0] = 1
    a[N][N] = 1
    
    b[0][0] = 0
    b[N][0] = 1
    

    for i in range(0, N-1):
        print(i)
        a[i+1][i] = 1
        a[i+1][i+1] = -2
        a[i+1][i+2] = 1
        
        
    a_inv = inv(a)
    
    result = np.matmul(a_inv, b)
    
    result_transpose = np.transpose(result)

    plt.figure(3)
    plt.title("Steady State Approach for Nonhomogeneous")

    plt.plot(list(range(len(result_transpose[0]))), result_transpose[0])
    plt.plot(result_transpose)

    
    '''
    -----------------------------------------------------------------------------------------
    part c: non-homogeneous transient (didnt end up working)
    
    '''
    
    # initialize temperature values at 0
    T = [T1] + [0 for i in range(N - 2)] + [T2]

    # run for 100000 steps
    for i in range(100000):
        
        new = [T1] + update(T, non_homo = False) + [T2]
        T = new
        
    plt.figure(4)
    plt.title("Transient Approach for Nonhomogeneous")
    
    # plot the final temperature distribution
    plt.plot(list(range(1,N+1)), T)
    
    
    
if __name__ == "__main__":
    main()


