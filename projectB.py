#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 13:27:35 2022

@author: ryanhuang
"""
import numpy as np
import matplotlib.pyplot as plt

'''
Converged scenarios with a decent amount of error
--------------------------------------------------------------------------------------------
For 5 Nodes, 200 Iterations, and without Under-Relaxing

Final Velocity: [1.3320748  1.7126676  2.39773464 3.99622439]
Final Pressure: [9.28135858 8.35710938 7.65768656 5.74913138 0.        ]

--------------------------------------------------------------------------------------------
For 5 Nodes, 200 Iterations, and Under-Relaxing

UR = 0.2, a = SU = SUL = SUR = 0.4, aL = 0.1, aR = 0.2, q = 0.1

Final Velocity: [1.35765631 1.74555811 2.44378135 4.07296892]
Final Pressure: [9.25340322 8.29324523 7.94834747 5.96707075 0.        ]

--------------------------------------------------------------------------------------------
For 50 nodes, 200 Iterations, and Under-Relaxing

UR = 0.15, a = 0.05, aL = 0.0001, aR = 0.02, SU = 0.99999999, SUL = 0.08, SUR = 0.3, q = 0.005

Velocity Error: 11.4975656% 

'''


# THESE PARAMETERS CAN BE CHANGED TO TEST DIFFERENT CASES
NUM_ITER = 200                                      # number of iterations
NODES = 50                                           # number of nodes to use
test_converge = False                                # break the loop if error reaches a low enough amount
under_relax = True                                   # must adjust the underrelaxation constants if True

# UNDER RELAXATION CONSTANTS
UR = 0.15                                            # for V and P 
a = 0.05                                             # for aP
aL = 0.0001                                            # for aP_right
aR = 0.02                                            # for aP_left
SU = 0.99999999                                            # for Su
SUL = 0.08                                           # for Su_left
SUR = 0.3                                           # for Su_right
q = 0.005                                            # for d


# THESE PARAMETERS ARE THE INITIAL CONDITIONS
POINTS = NODES - 1
DENS = 1                                            # density
P_inlet = 10                                        # inlet pressure
P_outlet = 0                                        # outlet pressure

# THESE VARIABLES REPRESENT THE INITIAL GUESSES MADE IN THE ALGORITHM
P_i = np.linspace(P_inlet, P_outlet, num=NODES)     # initial pressure guesses
area = np.linspace(0.5, 0.1, num=(NODES+POINTS))                                # area per point with linear spacing
AREA_nodes = np.array([area[i] for i in range(len(area)) if i % 2 == 0])        # area per node with linear spacing
AREA_points = np.array([area[i] for i in range(len(area)) if i % 2 !=0])        # area per point with linear spacing


# THESE VARIABLES ARE FOR CALCULATING PERCENT ERROR AT THE END 
exact_mass_flow_rate = 0.44721                      # exact mass flow rate
guessed_mass_flow_rate = 1                          # guessed mass flow_rate
exact_vel = np.array([exact_mass_flow_rate/(DENS * AREA_points[i]) for i in range(POINTS)])

     

def main():
    # variables for underrelaxation calculations
    UR_aP = 0
    UR_Su = 0
    UR_aP_left = 0
    UR_aP_right = 0
    UR_Su_right = 0
    UR_Su_left = 0
    UR_d = []
    
    # lists to track convergence
    aP_plt = []
    aP_R_plt = []
    aP_L_plt = []

    vel_error = []                                  # store the velocity error per iteration
    
    vel_mat = np.zeros((NODES-1,NODES-1))           # declare velocity matrix
    source_mat = np.zeros((NODES-1, 1))             # declare source matrix
    
    pressure_mat = np.zeros((NODES,NODES))          # declare pressure correction matrix
    temp_mat = np.zeros((NODES,1))                  # this matrix is used in conjunction with the pressure correction matrix
    
    
    V = np.array([guessed_mass_flow_rate/(DENS * AREA_points[i]) for i in range(POINTS)])

    print(f"Velocities: {exact_vel}\n")
    
    #---------------------------------------------------------------------
    #solve for the interior nodes
    
    for i in range(NUM_ITER):
        
        d = []                                          # parameter used in pressure correction
        
        for node in range(1, POINTS-1):                  # only consider the interior nodes and not the boundaries
            Fw = DENS * (0.5*(V[node-1] + V[node])) * AREA_nodes[node]
            Fe = DENS * (0.5*(V[node] + V[node+1])) * AREA_nodes[node+1]
            aW = Fw
            aE = 0
            
            # underrelaxation
            if i == 0 or under_relax == False:
                aP = aW + aE + Fe - Fw
                Su = (P_i[node] - P_i[node+1]) * AREA_points[node]
            else:
                aP = a*UR_aP + (1-a)*(aW + aE + Fe - Fw)
                Su = SU*UR_Su + (1-SU)*((P_i[node] - P_i[node+1]) * AREA_points[node])
            
            d.append(AREA_points[node]/aP)
            
            vel_mat[node][node-1] = -aW
            vel_mat[node][node] = aP
            source_mat[node] = Su
        
        #---------------------------------------------------------------------
        # solve for the boundary nodes
        #---------------------------------------------------------------------
        # leftmost boundary
        
        Fe_left = DENS * (0.5*(V[0] + V[1])) * AREA_nodes[1]
        
        uA = V[0] * (AREA_points[0]/AREA_nodes[0])
        Fw_left = DENS * uA * AREA_nodes[0]
        
        aW_left = aE_left = 0
        
        # underrelaxation
        if i == 0 or under_relax == False:
            aP_left = Fe_left + Fw_left * 0.5 * (AREA_points[0]/AREA_nodes[0])**2
            Su_left = (P_i[0] - P_i[1]) * AREA_points[0] + Fw_left * (AREA_points[0]/AREA_nodes[0]) * V[0]
        else: 
            aP_left = aL*UR_aP_left + (1-aL)*(Fe_left + Fw_left * 0.5 * (AREA_points[0]/AREA_nodes[0])**2)
            Su_left = SUL*UR_Su_left + (1-SUL)*((P_i[0] - P_i[1]) * AREA_points[0] + Fw_left * (AREA_points[0]/AREA_nodes[0]) * oldV[0])        
        
        
        vel_mat[0][0] = aP_left
        source_mat[0] = Su_left
        d.insert(0, AREA_points[0]/aP_left)
        
        #---------------------------------------------------------------------
        # rightmost boundary
        
        if i == 0:
            Fe_right = guessed_mass_flow_rate
        else:
            Fe_right = DENS * V[-1] * AREA_points[-1]
        
        Fw_right = DENS * (0.5*(V[-2] + V[-1])) * AREA_nodes[-2]
        aW_right = Fw_right
        aE_right = 0
        
        #underrelaxation
        if i == 0 or under_relax == False:
            aP_right = aW_right + aE_right + Fe_right - Fw_right
            Su_right = (P_i[-2] - P_i[-1]) * AREA_points[-1]
        else:
            aP_right = aR*UR_aP_right + (1-aR)*(aW_right + aE_right + Fe_right - Fw_right)
            Su_right = SUR*UR_Su_right + (1-SUR)*((P_i[-2] - P_i[-1]) * AREA_points[-1])
        
        vel_mat[-1][-2] = -aW_right
        vel_mat[-1][-1] = aP_right
        source_mat[-1] = Su_right
        d.append(AREA_points[-1]/aP_right)
        
        # underrelaxation for d
        temp_d = d

        if under_relax == True:
            if i != 0:
                for ele in range(len(d)):
                    d[ele] = q*UR_d[ele] + (1-q)*temp_d[ele]

        
        #---------------------------------------------------------------------
        #solve for the velocity field
        
        sol = np.linalg.solve(vel_mat, source_mat)
        
        #---------------------------------------------------------------------
        # pressure correction 
        
    
        for node in range(1, NODES-1):
            aW1 = DENS * d[node-1] * AREA_points[node-1]
            aE1 = DENS * d[node] * AREA_points[node]
            Fw1 = DENS * sol[node-1] * AREA_points[node-1]
            Fe1 = DENS * sol[node] * AREA_points[node]
            aP1 = aW1 + aE1 
            b = Fw1 - Fe1
            
            pressure_mat[node][node-1] = -aW1
            pressure_mat[node][node] = aP1
            pressure_mat[node][node+1] = -aE1
            
            
            temp_mat[node] = b
        
        # set the correction pressures at the boundaries equal to 0 by deleting the edges of the matrix

        #these two reassignments are two preserve the dimensions of pressure_mat and temp_mat 
        press_mat = pressure_mat
        temp1_mat = temp_mat
        
        press_mat = np.delete(press_mat, 0, axis = 0)
        press_mat = np.delete(press_mat, -1, axis = 0)
        press_mat = np.delete(press_mat, 0, axis = 1)
        press_mat = np.delete(press_mat, -1, axis = 1)
        temp1_mat = np.delete(temp1_mat, 0)
        temp1_mat = np.delete(temp1_mat, -1)
            
        
        press_corr = np.linalg.solve(press_mat, temp1_mat)
        #---------------------------------------------------------------------
        # correct pressures
        
        oldP = P_i

        for p in range(1,NODES-1):
            P_i[p] += press_corr[p-1]
        
        # add 0 to both ends of the pressure correction matrix to preserve length after the edges were deleted earlier
        press_corr = np.insert(press_corr, 0, 0)
        press_corr = np.append(press_corr, 0)
        
        #---------------------------------------------------------------------
        # correct velocities
        for v in range(0, POINTS):
            sol[v] += d[v] * (press_corr[v] - press_corr[v+1])
            
            
        #---------------------------------------------------------------------
        # solve for nodal pressure at A
        
        P_i[0] = P_inlet - 0.5 * DENS*((sol[0] * AREA_points[0])/AREA_nodes[0])**2

        #---------------------------------------------------------------------
        # underelaxation factor

        oldV = V
        
        
        for val in range(len(V)):
            V[val] = UR * oldV[val] + (1-UR) * sol[val]
            P_i[val] = UR * oldP[val] + (1-UR) * P_i[val]
            
        UR_aP = aP
        UR_Su = Su
        UR_aP_left = aP_left
        UR_aP_right = aP_right
        UR_Su_right = Su_right
        UR_Su_left = Su_left
        UR_d = d
        
        
        aP_plt.append(UR_aP)
        aP_R_plt.append(UR_aP_left)
        aP_L_plt.append(UR_aP_right)
        
        vel_error.append((V[0]/exact_vel[0] - 1)*100)
        
        
        # breaks the loop if an error reaches a certain point
        if test_converge:
            if len(vel_error) > 2:
                if abs(vel_error[-1]) < 0.01 or (abs(vel_error[-1] - vel_error[-2])) < 0.001:
                    print("converged?")
                    break

    # display the final velocities, pressures, and errors
    print(f"Final Velocity: {V}")
    print(f"Final Pressure: {P_i}\n")
    print(f"Final Error: {vel_error[-1]}%")
    
    
    # take the absolute value of all the errors to find the index of the minimum error
    min_err = [abs(err) for err in vel_error]
    min_err_index = min_err.index(min(min_err))
    
    print(f"Minimum Error: {vel_error[min_err_index]}")
    print(f"Iteration of Minimum Error : {min_err_index}")
    print(f"Last error difference: {vel_error[-1] - vel_error[-2]}")
    
    
    #plot the percent error per iteration
    plt.figure(0)
    plt.plot(list(range(len(vel_error))), vel_error)
    plt.xlabel("Iterations")
    plt.ylabel("Percent Error")
    plt.title(f"Velocity Percent Error at {NODES} nodes")
    
    
    #plot velocity distribution
    plt.figure(1)
    plt.title("Exact versus Numerical velocity distribution")
    plt.xlabel("Nodes")
    plt.ylabel("Velocity")
    plt.plot(list(range(len(exact_vel))), exact_vel, label="Exact")
    plt.plot(list(range(len(V))), V, label="Numerical")
    plt.legend()

    
if __name__ == "__main__":
    main()
