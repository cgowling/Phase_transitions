#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 14:04:14 2018

@author: cg411
Looking at the fluid shell function in the bubble module, 
trying to recreate figure 5 in the energy budget paper .
"""

import bubble as b
import matplotlib.pyplot as plt
import numpy as np 


#%%

v_w_input_arr = np.linspace(0.05,0.8,7)



#getting plot for deflagration to detontation V,xi. Also enthalpy plots

fig, axes = plt.subplots(2, 1)
ax1 = plt.subplot2grid((1, 2), (0, 0))
ax2 = plt.subplot2grid((1, 2), (0, 1))

aN0 = 0.3
for i in  np.arange(0,len(v_w_input_arr)):
        
    a_plus = b.find_alpha_plus(v_w_input_arr[i],aN0)
    v_f, w_f,xiM = b.fluid_shell(v_w_input_arr[i],a_plus)
    
    ax1.plot(xiM,v_f)
    ax2.plot(xiM,w_f)
  
#ax1.set_xlabel('xi')
#ax1.set_ylabel('v(xi)')      
#
#  
#ax2.set_xlabel('xi')
#ax2.set_ylabel('w(xi)')    
##ax2.
#plt.autoscale
#
#
#ax1.plt.sca(axes[0])


plt.show
#
##%% 
##getting fig6
#aN = 0.1
##xi_w = np.linspace(0,0.9,10)
#xi = b.make_xi_array()
#
#
#v_w_input_arr = np.linspace(0.05,0.8,7)
#
#
#
#
##%%
##getting plot Fig.8 energy coefficient    
#aN = [0.01,0.03,0.1,0.3,1,3,10]
##xi = b.make_xi_array()
#
#
#v_w_input_arr = np.linspace(0.05,0.8,30)
#kappa = np.zeros((7,30))
#for i in np.arange(0,len(aN)):
#    
#    kappa[i,:] = b.get_kappa_arr(v_w_input_arr,aN[i])
#
#xi_w = np.linspace(0,0.8,30)
#plt.plot(xi_w,kappa[0:4,:])
#for i in  np.arange(0,len(v_w_input_arr)):
#        
#    a_plus = b.find_alpha_plus(v_w_input_arr[i],aN0)
#    v_f, w_f,xiM = b.fluid_shell(v_w_input_arr[i],a_plus)
##kappa = np.zeros(len(v_w_input_arr),len(v_w_input_arr))    
#for i in  np.arange(0,len(v_w_input_arr)):
#        
#  kappa = b.get_kappa_arr(v_w_input_arr,aN[i])
#  plt.plot(xiM, kappa)
#  print(b.get_kappa_arr(v_w_input_arr,aN[i]))
#                                 
#    plt.plot(xiM,kappa[i])
#    plt.plot(xiM,v_f)
#    
#    plt.xlabel('xi')
#   
##    ax2.plot(xiM,w_f)
#
#    plt.ylabel('v(xi)')
#
#plt.show()
#Value of v_wall how to select?? 