# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 14:24:17 2022

@author: Administrator
"""

#!/usr/bin/env python
# Eclipse SUMO, Simulation of Urban MObility; see https://eclipse.org/sumo
# Copyright (C) 2009-2017 German Aerospace Center (DLR) and others.
# This program and the accompanying materials
# are made available under the terms of the Eclipse Public License v2.0
# which accompanies this distribution, and is available at
# http://www.eclipse.org/legal/epl-v20.html

# @file    runner.py
# @author  Lena Kalleske
# @author  Daniel Krajzewicz
# @author  Michael Behrisch
# @author  Jakob Erdmann
# @date    2009-03-26
# @version $Id$

from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import optparse
import subprocess
import random
import math
import matplotlib.pyplot as plt
from scipy.linalg import expm
from numpy.linalg import inv
import traci
import numpy as np 
import numpy.matlib
from scipy.io import loadmat
from numpy import linalg as LA
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import csv

# we need to import python modules from the $SUMO_HOME/tools directory
os.environ['SUMO_HOME'] = "/opt/homebrew/opt/sumo/share/sumo"

State_pre = np.zeros(11)
State_now = np.zeros(11)
Tra_Ld = np.array([[0,0]])
Time = np.array(0.01)
cross_Ld = 0 # record the times the AV crosses the 0 degree axis
cross_Fl = 0

try:
    sys.path.append(os.path.join(os.path.dirname(
        __file__), '..', '..', '..', '..', "tools"))  # tutorial in tests
    sys.path.append(os.path.join(os.environ.get("SUMO_HOME", os.path.join(
        os.path.dirname(__file__), "..", "..", "..")), "tools"))  # tutorial in docs
    from sumolib import checkBinary  # noqa
except ImportError:
    sys.exit(
        "please declare environment variable 'SUMO_HOME' as the root directory of your sumo installation (it should contain folders 'bin', 'tools' and 'docs')")


def generate_routefile():
     # demand per second from different directions
     with open("CircularRoad.rou.xml", "w") as routes: 
         # laneChangeModel= "SL2015"
         print("""<routes>
         <vType id="auto" accel="500" decel="500" length="5" minGap="1" maxSpeed="30" guiShape="passenger" color="1,0,0"/>
         <vType id="human" accel="5" decel="4.5"  length="5" minGap="1" maxSpeed="30" guiShape="passenger" color="0,1,0"/>
         <vType id="human1" accel="5" decel="4.5"  length="5" minGap="1" maxSpeed="30" guiShape="passenger" color="0,1,0"/>
         <vType id="leader" accel="5" decel="4.5"  length="5" minGap="1" maxSpeed="30" guiShape="passenger" color="0,1,0"/>
         <route id="right" edges="edge1 edge2"/>""", file=routes)

         print('<vehicle id="Ld" type="human" route="right" depart="0" departSpeed="20" departLane="1" departPos="50" />', file=routes)
         print('<vehicle id="Fl" type="auto" route="right" depart="0" departSpeed="20" departLane="0" departPos="20" />', file=routes)        
         
         print("</routes>", file=routes)

def get_ud(X): 
    # Parameters of the following vehicle
    M = 1360.0
    Iz = 1993.0
    lf = 1.45
    lr = 1.06
    ls = 0.71
    Ca = 0.5
    k0 = 460.0/0.33
    Cf = 1.51*10**5
    Cr = 1.46*10**5


    # Get the state of the leading vehicles
    x_L =  X[0]
    y_L =  X[1]
    phi_L =  X[2] 
    V_L =  X[3] 
    omega_L =  X[4]
    
    
    # Compute the solution to the regulator equation 
    g_v = np.array([[2*k0/M,       0,         2*k0/M],
                       [0,            Cf/M,      0],
                       [- ls*2*k0/Iz, Cf*lf/Iz,  ls*2*k0/Iz]])
    
    ud = inv(g_v)@np.array([[Ca/M*V_L**2],
                            [V_L*omega_L - (Cr*lr - Cf*lf)/M*omega_L/V_L],
                            [(Cf*lf**2 + Cr*lr**2)/Iz*omega_L/V_L]])

    return ud

def multichoose(k, objects):
    """n multichoose k multisets from the list of objects.  n is the size of
    the objects."""
    j,j_1,q = k,k,k  # init here for scoping
    r = len(objects) - 1
    a = [0 for i in range(k)] # initial multiset indexes
    while True:
        yield [objects[a[i]] for i in range(0,k)]  # emit result
        j = k - 1
        while j >= 0 and a[j] == r: j -= 1
        if j < 0: break  # check for end condition
        j_1 = j
        while j_1 <= k - 1:
            a[j_1] = a[j_1] + 1 # increment
            q = j_1
            while q < k - 1:
                a[q+1] = a[q] # shift left
                q += 1
            q += 1
            j_1 = q

def get_error(X,kappa):

    
    
    # Get the state of the vehicles
    x_L = X[0]
    y_L = X[1]
    phi_L = X[2] 
    V_L = X[3]
    omega_L = X[4]
    
    x = X[5]
    y = X[6]
    phi = X[7]
    V_x = X[8]
    V_y = X[9]
    omega = X[10]
    
    # parameters of the path
    ts = 1
    ds = 5
    d = ts*V_x + ds
    gamma = math.atan(kappa*d)
    
    if kappa == 0:
        s = 0
    else:
        s = (-1 + math.sqrt(1+kappa**2*d**2))/kappa
        
    # Get the error state
    beta = math.atan(V_y/V_x)
    e1 = x_L + s*math.sin(phi_L) - x - d*math.cos(phi)
    e2 = y_L - s*math.cos(phi_L) - y - d*math.sin(phi)
    e3 = phi_L - (phi + gamma)
    e4 = V_L - V_x
    e5 = -V_y
    e6 = omega_L - omega
    z1 =  math.cos(phi+gamma)*e1 + math.sin(phi+gamma)*e2;
    z2 = -math.sin(phi+gamma)*e1 + math.cos(phi+gamma)*e2;

    return [z1,z2,e3,e4,e5,e6]

def get_basis_ue(X, kappa):
    
    V_x = X[8]
    e = get_error(X,kappa)
    basis_ue = np.array([])
    
    cs1 = np.array([1-math.cos(e[2]), math.sin(e[2])])
    
    poly2 = np.array([])   
    for comb in multichoose(2, e):
        poly2 = np.append(poly2, np.prod(comb))
    basis_ue = np.append(basis_ue,poly2)
    
    poly1 = np.array([])  
    for comb in multichoose(1, e):
        poly1 = np.append(poly1, np.prod(comb))
    basis_ue = np.append(basis_ue,poly1)
    

    basis_ue = np.append(basis_ue, np.kron(poly1,cs1))
    # basis_ue = np.append(basis_ue, np.kron(poly1,V_x))
    # basis_ue = np.append(basis_ue, np.kron(cs1,V_x))
    basis_ue = np.append(basis_ue, np.array([e[3]/V_x, e[4]/V_x, e[5]/V_x]))
    
    return basis_ue
    
def get_ue0(X, kappa):
    #  Parameters of the following vehicle
    M = 1360.0
    Iz = 1993.0
    lf = 1.45
    lr = 1.06
    ls = 0.71
    Ca = 0.5
    k0 = 460.0/0.33
    Cf = 1.51*10**5
    Cr = 1.46*10**5
    
    K1 = 1*np.identity(3)
    K2 = 5*np.identity(3)
        
    
    # Get the state of the vehicles
    x_L =  X[0]
    y_L =  X[1]
    phi_L =  X[2] 
    V_L =  X[3]
    omega_L =  X[4]
    
    x =  X[5]
    y =  X[6]
    phi =  X[7]
    V_x =  X[8]
    V_y =  X[9]
    omega =  X[10]
    
    # parameters of the path
    ts = 1
    ds = 5
    d = ts*V_x + ds
    gamma = math.atan(kappa*d)
    
    if kappa == 0:
        s = 0
    else:
        s = (-1.0 + math.sqrt(1+kappa**2*d**2))/kappa
    
    # calculate error states
    e1 = x_L + s*math.sin(phi_L) - x - d*math.cos(phi)
    
    # Get the error state
    e2 = y_L - s*math.cos(phi_L) - y - d*math.sin(phi)
    e3 = phi_L - (phi + gamma)
    e4 = V_L - V_x
    e5 = -V_y
    e6 = omega_L - omega
    
    e123 = np.array([[e1],[e2],[e3]])
    e456 = np.array([[e4],[e5],[e6]])
    
    
    # Calculate Controller
    f_ep = np.array([[V_L*math.cos(phi_L) + s*math.cos(phi_L)*omega_L - V_L*math.cos(phi_L - gamma - e3) + d*math.sin(phi_L - gamma - e3)*omega_L],
                         [V_L*math.sin(phi_L) + s*math.sin(phi_L)*omega_L - V_L*math.sin(phi_L - gamma - e3) - d*math.cos(phi_L - gamma - e3)*omega_L],
                         [0]])   
    f_ep_dot = np.array([[-V_L*math.sin(phi_L)*omega_L - s*math.sin(phi_L)*omega_L**2 + V_L*math.sin(phi_L-gamma-e3)*(omega_L - e6) + d*math.cos(phi_L-gamma-e3)*omega_L*(omega_L - e6)],
                 [V_L*math.cos(phi_L)*omega_L + s*math.cos(phi_L)*omega_L**2 - V_L*math.cos(phi_L-gamma-e3)*(omega_L - e6) + d*math.sin(phi_L-gamma-e3)*omega_L*(omega_L - e6)],
                 [0]])
    g_ep = np.array([[math.cos(phi_L - gamma - e3), -math.sin(phi_L - gamma - e3), -d*math.sin(phi_L - gamma - e3)],
                     [math.sin(phi_L - gamma - e3),  math.cos(phi_L - gamma - e3),  d*math.cos(phi_L - gamma - e3)],
                     [0,                       0,                       1]])
    g_ep_inv = np.array([[math.cos(phi_L-gamma-e3), math.sin(phi_L-gamma-e3),  0],
                         [-math.sin(phi_L-gamma-e3), math.cos(phi_L-gamma-e3), -d],
                         [0,                       0,                       1]])
    g_ep_inv_dot = np.array([[-math.sin(phi_L-gamma-e3),  math.cos(phi_L-gamma-e3),   0],
                             [-math.cos(phi_L-gamma-e3),  -math.sin(phi_L-gamma-e3),  0],
                             [0,           0,          0]])*(omega_L - e6)   
    
    e123_dot = f_ep + g_ep@e456    
    alphae = g_ep_inv@(-f_ep - K1@e123)
    alphae_dot = g_ep_inv_dot@(-f_ep - K1@e123) + g_ep_inv@( -f_ep_dot - K1@e123_dot)
    
    f_ev = np.array([[-V_y*omega + Ca/M*V_x**2 - Ca/M*V_L**2],
                     [V_x*omega - V_L*omega_L + (Cr + Cf)/M*V_y/V_x + (Cr*lr - Cf*lf)/M*(omega_L/V_L - omega/V_x)],
                     [-(Cr*lr - Cf*lf)/Iz*V_y/V_x + (Cf*lf**2 + Cr*lr**2)/Iz*(omega/V_x- omega_L/V_L)]])
                 
    g_v = np.array([[2*k0/M,       0,         2*k0/M],
                    [0,            Cf/M,      0],
                    [- ls*2*k0/Iz, Cf*lf/Iz,  ls*2*k0/Iz]])
        
    ue = inv(g_v)@(f_ev - alphae_dot + K2@(e456 - alphae) + np.transpose(g_ep)@e123)
    
    return ue


def model(X,t,kappa):
# DINAMICS the dynamic model of the leading vehicle and the following vehicle
# X[0] = x_L, X[1] = y_L, X[2] = \phi_L, X[3] = V_L, X[4] = \omega_L
# X[5] = x, X[6] = y, X[7] = \phi, X[8] = V_x, X[9] = V_y, X[10] = omega
  
    # Parameters of the following vehicle
    M = 1360.0
    Iz = 1993.0
    lf = 1.45
    lr = 1.06
    ls = 0.71
    Ca = 0.5
    k0 = 460.0/0.33
    Cf = 1.51*10**5
    Cr = 1.46*10**5
    
    # Get the state of the vehicles
    x_L =  X[0]
    y_L =  X[1]
    phi_L =  X[2]
    V_L =  X[3]
    omega_L =  X[4]

    x =  X[5]
    y =  X[6]
    phi =  X[7] 
    V_x =  X[8] 
    V_y =  X[9]
    omega =  X[10]

    # Kinematics of the leading vehicle
    x_L_dot = V_L*math.cos(phi_L)
    y_L_dot = V_L*math.sin(phi_L)
    phi_L_dot = omega_L
    V_L_dot = 0
    omega_L_dot = 0
    
    # get the control input
    ud = get_ud(X)
    if t<0:
        ue = get_ue0(X,kappa)
    else:
        Su = loadmat('Su.mat')
        Su = Su['Su']
        Su = np.matrix(Su)
        basis_ue = get_basis_ue(X, kappa)
        basis_ue = np.reshape(basis_ue,(len(basis_ue),1))
        ue = Su@basis_ue
    
    u = ud + ue
    u1 = u[0][0]
    u2 = u[1][0]
    u3 = u[2][0]
    
    # Dynamics of the following vehicle    
    x_dot = V_x*math.cos(phi) - V_y*math.sin(phi)
    y_dot = V_x*math.sin(phi) + V_y*math.cos(phi)
    phi_dot = omega
    V_x_dot = V_y*omega - Ca/M*V_x**2 + 2*k0/M*u1 + 2*k0/M*u3
    V_y_dot = -V_x*omega - (Cf + Cr)/M*V_y/V_x + (Cr*lr - Cf*lf)/M*omega/V_x  + Cf/M*u2
    omega_dot = (Cr*lr - Cf*lf)/Iz*V_y/V_x - (Cf*lf**2 + Cr*lr**2)/Iz*omega/V_x - ls*2*k0/Iz*u1 + ls*2*k0/Iz*u3 + Cf*lf/Iz*u2    
    
    dX = [x_L_dot, y_L_dot, phi_L_dot, V_L_dot, omega_L_dot, x_dot, y_dot, phi_dot, V_x_dot, V_y_dot, omega_dot]
    
    return dX


def run():
     """execute the TraCI control loop"""
     dt = 0.01
     step = 0
     N1 = 3000
     f = open('./State_tra.csv', 'w')
     writer = csv.writer(f)
     global Time, State_pre, cross_Ld, cross_Fl, Tra_Ld, State_now
     
     while step < N1:
         step += 1
         traci.simulationStep()
         traci.vehicle.setSpeedMode('Ld', 0)
         traci.vehicle.setSpeedMode('Fl', 0)
         
         traci.vehicle.setLaneChangeMode('Ld', 1)
         traci.vehicle.setLaneChangeMode('Fl', 1)  
         
         Time = np.append(Time,traci.simulation.getTime())
         
         # Measure the state of the vehicles
         Pos_Ld_y = traci.vehicle.getPosition('Ld')[1]
         Pos_Ld_x = traci.vehicle.getPosition('Ld')[0]
         Tra_Ld = np.append(Tra_Ld,[[Pos_Ld_x,Pos_Ld_y]],axis=0)
         
         Pos_Fl_y = traci.vehicle.getPosition('Fl')[1]
         Pos_Fl_x = traci.vehicle.getPosition('Fl')[0]
    
         
         kappa = 1/(51.6) #curvature of the road
         if Time[-1]>1:
             dx_dt = np.gradient(Tra_Ld[-100:-1, 0])
             dy_dt = np.gradient(Tra_Ld[-100:-1, 1])
             ds_dt = np.sqrt(dx_dt * dx_dt + dy_dt * dy_dt)
             d2s_dt2 = np.gradient(ds_dt)
             d2x_dt2 = np.gradient(dx_dt)
             d2y_dt2 = np.gradient(dy_dt)
             curvature = np.abs(d2x_dt2 * dy_dt - dx_dt * d2y_dt2) / (dx_dt * dx_dt + dy_dt * dy_dt)**(1.5)
             kappa0 = np.mean(curvature)
             print(1/kappa0)
         
         vel_Ld_x = traci.vehicle.getSpeed('Ld')
         
         angle_Ld = traci.vehicle.getAngle('Ld')
         angle_Ld = -math.radians(angle_Ld) + math.pi/2 + 2*math.pi*cross_Ld
         if abs(angle_Ld-State_pre[2])>math.pi*2-0.5: # when the yaw angle cross the 2*pi, it will be set to zero by sumo, leading to the discontinuity of yaw angle
             cross_Ld += 1 
             angle_Ld += 2*math.pi
         #omg_Ld = (angle_Ld - State_pre[2])/(Time[-1]-Time[-2]) if (Time[-1]-Time[-2])>0 else State_pre[4]
         omg_Ld = vel_Ld_x*kappa    
         
         
         vel_Fl_x = traci.vehicle.getSpeed('Fl')
         vel_Fl_y = traci.vehicle.getLateralSpeed('Fl')        

         
         angle_Fl = traci.vehicle.getAngle('Fl')
         angle_Fl = -math.radians(angle_Fl) + math.pi/2 + 2*math.pi*cross_Fl
         if abs(angle_Fl-State_pre[7])>math.pi*2-0.5: # when the yaw angle cross the 2*pi, it will be set to zero by sumo, leading to the discontinuity of yaw angle
             cross_Fl += 1 
             angle_Fl += 2*math.pi
         omg_Fl = (angle_Fl - State_pre[7])/(Time[-1]-Time[-2]) if (Time[-1]-Time[-2])>0 else State_pre[10]   
         #omg_Fl = math.sqrt(vel_Fl_x**2 + vel_Fl_y**2)*kappa

         
         if step<=1:
             X0 = [Pos_Ld_x, Pos_Ld_y, angle_Ld, vel_Ld_x, omg_Ld, Pos_Fl_x, Pos_Fl_y, angle_Fl, vel_Fl_x, vel_Fl_y, omg_Fl]
         else:
             X0 = [Pos_Ld_x, Pos_Ld_y, angle_Ld, vel_Ld_x, omg_Ld, State_now[5], State_now[6], State_now[7], State_now[8], State_now[9], State_now[10]]
             # X0 = [Pos_Ld_x, Pos_Ld_y, angle_Ld, vel_Ld_x, omg_Ld, Pos_Fl_x, Pos_Fl_y, angle_Fl, vel_Fl_x, State_now[9], State_now[10]]
             # X0 = [State_now[0], State_now[1], State_now[2], State_now[3], State_now[4], State_now[5], State_now[6], State_now[7], State_now[8], State_now[9], State_now[10]]
         # error = get_error(X0, kappa)
         # print(error)
         
         #Compute the dynamics and kinematics of the vehicles
         State_pre = np.array(X0)
         t = np.linspace(Time[-1],Time[-1]+dt,11)
         X = odeint(model,X0,t,atol=1e-10, rtol=1e-10,args=(kappa,))
         error = get_error(X[-1,:], kappa)
         State_now = np.array(X[-1,:])
         
         #Set the position and velocity of the vehicles
         omega_L = 20*kappa
         theta = omega_L*Time[-1]+math.pi/4
         x_L = 50.0 + 1/kappa*math.cos(theta)
         y_L = 50.0 + 1/kappa*math.sin(theta)
         traci.vehicle.moveToXY(vehID='Ld',edgeID="", lane=-1, x=x_L, y=y_L, angle=-math.degrees(theta), keepRoute=2)
         traci.vehicle.setSpeed('Ld', 20)
         
         traci.vehicle.moveToXY(vehID='Fl',edgeID="", lane=-1, x=X[-1,5], y=X[-1,6], angle=-math.degrees(X[-1,7])+90.0, keepRoute=2)
         # traci.vehicle.setSpeed('Fl', X[-1,8])

        
         writer.writerow(State_pre)
     f.close()
         
def get_options():
     optParser = optparse.OptionParser()
     optParser.add_option("--nogui", action="store_true",
                          default=False, help="run the commandline version of sumo")
     options, args = optParser.parse_args()
     return options


# this is the main entry point of this script
if __name__ == "__main__":
     options = get_options()

     # this script has been called from the command line. It will start sumo as a
     # server, then connect and run
     if options.nogui:
         sumoBinary = checkBinary('sumo')
     else:
         sumoBinary = checkBinary('sumo-gui')

     # first, generate the route file for this simulation
     generate_routefile()

     # this is the normal way of using traci. sumo is started as a
     # subprocess and then the python script connects and runs
     traci.start([sumoBinary, "-c", "CircularRoad.sumocfg"])
     run()