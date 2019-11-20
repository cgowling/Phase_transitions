#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 09:41:20 2019

@author: cg411
"""

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl

import pttools.gw_code.ssmtools as s
import pttools.bubble.bubble as b

import GW_tools as g 

from scipy.optimize import curve_fit
    
mpl.rcParams.update({"font.size": 30})



            
 
    
def fit_gw_ssm(z,params,model,lower_bounds,upper_bounds, initial_guess):
    
    
    om_gw_ssmtools = s.power_gw_scaled(z, params, npt = [5000,1000])
    
    
    initial_guess[1] = max(om_gw_ssmtools)
    sigma = max(om_gw_ssmtools)* np.ones(om_gw_ssmtools.shape)
    param_estimate , _ = curve_fit(model,z,om_gw_ssmtools, bounds =(lower_bounds, upper_bounds),p0 = initial_guess, sigma = sigma )
    return om_gw_ssmtools, param_estimate
    

     
def compare_fit(z, om_ssm,vw_indices ,alphas,PEs): 
    
    for j , vw_idx in enumerate(vw_indices):  
        plt.figure()
        
        for i,alpha in enumerate(alphas):
            plt.loglog(z,om_ssm[i,vw_idx ,:], label= alpha)
            
#            zb,Om_peak,rb
            om_fit = g.gw_spec_ssm_rescale(z,PEs[i,vw_idx,0],PEs[i,vw_idx,1],PEs[i,vw_idx,2] )
            plt.loglog(z,om_fit,'--')
        plt.title(r'$v_w= {}$'.format(vws[vw_idx]))
        get_plot_setup_om_z()


    

def polyfitter(vws,alphas,PEs ,deg, n):
    
        
            
    p_out = np.zeros([len(alphas), deg+1]) 
    deg_int = np.linspace(0,deg,deg +1) # NOTE polyfit returns in descending order !!!!!
    plt.figure()
    for i,alpha in enumerate(alphas):
        p_out[i,:] = np.polyfit(vws,PEs[i,:,n],deg)
        p = np.poly1d(p_out[i,:])
        poly_fit_rp = p(vws)
    
        
        if n == 1 :
            plt.semilogy(vws,poly_fit_rp,'--')
            plt.semilogy(vws,PEs[i,:,n],label = alpha)            
        else :
            plt.plot(vws,poly_fit_rp,'--')
            plt.plot(vws,PEs[i,:,n],label = alpha)
        
    plt.legend(title = 'alpha') 
     
    plt.figure()  
    for j, d_int in enumerate(deg_int):
            plt.plot(alphas,p_out[:,j],label = deg -d_int)
    
    
    
    plt.legend(title = 'degree of polynomial coef')    
    
    

def get_plot_setup_om_z():   
        plt.legend()
        plt.xlabel(r'$z$')
        plt.ylabel(r'$\Omega$')

        plt.legend(title= r'$\alpha$')
        plt.tight_layout()
             
    

def get_plot_setup():     
        
  
    plt.xlabel(r'$v_w$')
    plt.legend(title= r'$\alpha$')
    plt.grid(True, which='both',alpha=0.5)
    plt.xlim( (0.4,0.9))
    plt.tight_layout()
    return plt.gca()


def plot_parameter_estimates(vws,alphas,param_estimate,y_label):
    plt.figure()
    for i, alpha in enumerate(alphas): 
        plt.plot(vws,param_estimate[i], label=  alpha )
        
    plt.ylabel(y_label) 
    get_plot_setup()
    
def plot_parameter_estimates_semilogy(vws,alphas,param_estimate,y_label):
    plt.figure()
    for i, alpha in enumerate(alphas): 
        plt.semilogy(vws,param_estimate[i], label=  alpha )
        
    plt.ylabel(y_label) 
    get_plot_setup()


def alpha_vw_to_rb(alpha, vw):
    a0 =  -32.892* alpha + 6.18 #11.4
    a1 = 173.9*alpha -26.7
    a2 = -287.75*alpha +	37.4
    a3 = 150.3* alpha -16.58
    
    return a0 + a1 * vw + a2 * vw**2 + a3 * vw**3



def alpha_vw_to_omp(alpha, vw):
    a0 = 0.0209015 * alpha**2  -0.00147626	* alpha +1.70934e-05 
    a1 = -0.0401596* alpha**2 +0.00365655*alpha -4.42731e-05
    a2 = 0.00935957*alpha**2 -0.002548	*alpha +	3.39144e-05
    a3 =0.0117989 *alpha**2 +0.000273503*alpha -5.74102e-06
    


    
    return a0 + a1 * vw + a2 * vw**2   + a3* vw**3



def omp_vary_vw(alpha, vws):
    omp_fit = np.zeros(len(vws))
    print('called')
    for i, v in enumerate(vws):
        omp_fit[i] = alpha_vw_to_omp(alpha, v)
        
    return omp_fit
        
#%%
    
# Iterating over alpha and vw parameter space performing curve fit 

z = np.logspace(-1,2, num = 1000 , base = 10) #needs to be long enough so accurate within ssmtools 

vws = np.linspace(0.4,0.9,100) # below 0.4 not many (if any simulations done) 
alphas = [0.01,0.02,0.04,0.06, 0.07,0.1]#np.linspace(0.01,0.1)#

n_vw = len(vws)
n_alpha = len(alphas)


param_estimates = np.zeros([n_alpha,n_vw, 3])
om_gw_ssm = np.zeros([n_alpha,n_vw,len(z)])


# Bounds for parameterisation #4 rescaled 
#z,zb,Om_peak,rb 
lower_b = [0,0,0]
upper_b = [5,np.inf,1]

initial_guess = [0.5,0,0.5]





#%%

for i,alpha in enumerate(alphas):
    print(alpha)
    
    for j,vw in enumerate(vws): 
        
        params = [vw,alpha,'exponential',(1,)]
        
        om_gw_ssm[i,j],param_estimates[i,j,:] = fit_gw_ssm(z,params,g.gw_spec_ssm_zb,lower_b,upper_b,initial_guess)
    

#  using the z for internal integration , ideally would want this at a higher res, for now just use more zs ( could have them more densely packed around cs)
#switch to a diffferent array around the speed of sound 
#or add more values to the array as you loop 


#np.save('PEs',param_estimates)

#%%
vw_indices = [0,10,20,30,40,50,60,70,80,90]


compare_fit(z, om_gw_ssm,vw_indices ,alphas,param_estimates)
        
#%%

param_estimates   = np.load ( 'PEs.npy' ) 
om_tilde = np.load('om_tilde.npy')
#%%
#dependency on the break ratio rb


plot_parameter_estimates(vws,alphas,param_estimates[:,:,2],r'$r_b$')

    
plot_parameter_estimates(vws,alphas,param_estimates[:,:,0],r'$z_b$')   

plot_parameter_estimates(vws,alphas,param_estimates[:,:,0]/param_estimates[:,:,2],r'$z_p$')
    

plot_parameter_estimates_semilogy(vws,alphas,param_estimates[:,:,1],r'$\Omega_p$')
    
plot_parameter_estimates(vws,alphas,param_estimates[:,:,1],r'$\Omega_p$')       


#%%

#Approximate dependence on om_gw tilde - HnR* ???


om_tilde = np.zeros([n_alpha, n_vw])

plt.figure()
for i,alpha in enumerate(alphas):

    
    for j,vw in enumerate(vws): 
#                                           omp,vw,alpha,rb
        om_tilde[i,j] = g.om_to_om_tilde(param_estimates[i,j,1],vw,alpha, param_estimates[i,j,2])
    plt.semilogy(vws,om_tilde[i,:], label=  alpha )
        
plt.ylabel(r'$\tilde{\Omega}_{gw}$')
get_plot_setup()





#%% 
# finding analytic approximation for dependencies on vw and alpha
deg = 3
coefs = np.zeros([len(alphas), deg+1]) 
deg_int = np.linspace(0,deg,deg +1)

plt.figure()
for i,alpha in enumerate(alphas):

    coefs[i,:] = np.polyfit(vws,param_estimates[i,:,2],deg)
    p = np.poly1d(coefs[i,:])
    poly_fit_rb = p(vws)

    
    
    plt.plot(vws,poly_fit_rb,'--')
    plt.plot(vws,param_estimates[i,:,2],label = alpha)
    
plt.legend(title = 'alpha') 
 
plt.figure()  
for j, d_int in enumerate(deg_int):
        plt.plot(alphas,coefs[:,j],label = deg- d_int)



plt.legend(title = 'degree of polynomial coef')  

#%%

plt.figure()
coefs_alpha_dep = np.zeros( [deg+1,2])


for i in range(0,4):
    print(i)
    idx = deg-i

    coefs_alpha_dep[i,:] =np.polyfit(alphas,coefs[:,idx],1)

    p1 = np.poly1d(coefs_alpha_dep[i,:])
    poly_fit_alpha = p1(alphas)

    plt.plot(alphas,coefs[:,i],label = i)
    plt.plot(alphas, poly_fit_alpha, '--')
plt.legend()

#%%
deg = 3
coefs_omp = np.zeros([len(alphas), deg+1]) 



plt.figure()
for i,alpha in enumerate(alphas):

    coefs_omp[i,:] = np.polyfit(vws,param_estimates[i,:,1],deg)
    p_omp = np.poly1d(coefs_omp[i,:])
    poly_fit_omp = p_omp(vws)

    
    
    plt.semilogy(vws,poly_fit_omp,'--')
    plt.semilogy(vws,param_estimates[i,:,1],label = alpha)
    
plt.legend(title = 'alpha') 
 
plt.figure()  
for j, d_int in enumerate(deg_int):
        plt.plot(alphas,coefs_omp[:,j],label = deg- d_int)



plt.legend(title = 'degree of polynomial coef omp')  
#%%

plt.figure()
deg2 = 2
coefs_alpha_dep_omp = np.zeros([4,3]) #[deg2+1,3])

for i in range(0,4):
       idx = deg-i 
   
       coefs_alpha_dep_omp[i,:] =np.polyfit(alphas,coefs_omp[:,idx],deg2)
       p1_omp = np.poly1d(coefs_alpha_dep_omp[i,:])
       poly_fit_alpha_omp = p1_omp(alphas) 
       p1_omp = np.poly1d(coefs_alpha_dep_omp[i,:])
       poly_fit_alpha_omp = p1_omp(alphas)   
        
       plt.plot(alphas,coefs_omp[:,i],label = i)
       plt.plot(alphas, poly_fit_alpha_omp,'--')




#%%
#Comparing fit and ssmtools
plt.figure()

vw_tests = [0.4,0.5,0.6,0.7,0.8,0.9]
#alpha = alphas 


zb = 1 
for i,v in enumerate(vw_tests):
    plt.figure()
    for j,a in  enumerate(alphas):

        omp = alpha_vw_to_omp(alphas[j], vw_tests[i])
        
        rb = alpha_vw_to_rb(alphas[j], vw_tests[i])
        
        om_ssm_obs = -g.gw_spec_ssm_rescale(z,zb,-omp ,rb )
        
        
        params = [vw_tests[i],alphas[j],'exponential',(1,)]
        
        om_ssmtools = s.power_gw_scaled(z, params, npt = [5000,1000])
        
        
            
        plt.loglog(z,om_ssm_obs, '--' ,label = (vw_tests,  alphas))
            
        plt.loglog(z,om_ssmtools )


plt.legend()
#get_plot_setup_om_z()

# %%
#plt.figure()
#    
#z_test = np.polynomial.polynomial.Polynomial.fit(vws[40:99],param_estimates[0,40:99,2],7 )
#coefs = z_test.coef
#p = np.poly1d(coefs)
#poly_fit_rp = p(vws[40:99])
#plt.plot(vws[40:99],poly_fit_rp,'--')
#plt.plot(vws[40:99],param_estimates[0,40:99,2])




#
#def compare_fit_2p(z, om_ssm,vw_indices ,alphas,PEs): 
#    
#    for j , vw_idx in enumerate(vw_indices):  
#        plt.figure()
#        
#        for i,alpha in enumerate(alphas):
#            plt.loglog(z,om_ssm[i,vw_idx ,:], label= alpha)
#    
#            om_fit = g.gw_spec_ssm_const_zb(z,PEs[i,vw_idx,0],PEs[i,vw_idx,1] )
#            plt.loglog(z,om_fit,'--')
#    get_plot_setup_om_z()

#def fit_gw_ssm_2p(z,params,model,lower_bounds,upper_bounds, initial_guess):
#    
#    
#    om_gw_ssm = s.power_gw_scaled(z, params, npt = [5000,1000])
#    
#    
#    initial_guess[0] = max(om_gw_ssm)
##    print(initial_guess)
#    sigma = max(om_gw_ssm)* np.ones(om_gw_ssm.shape)
#    param_estimate , _ = curve_fit(model,z,om_gw_ssm, bounds =(lower_bounds, upper_bounds),p0 = initial_guess, sigma = sigma )
#    return om_gw_ssm, param_estimate
#

#
#
#def fit_gw_ssm_1(z,params,model):
#    
#        om_gw_ssm = s.power_gw_scaled(z, params, npt = [5000,1000])
#        
#        # curve_fit set up 
#        
#        lower_bounds = [0,0,0]
#        upper_bounds = [np.inf,np.inf,1] # change back if doing rp !!!
#        
#        
#        initial_guess = [1,max(om_gw_ssm),0.5]
#        sigma = max(om_gw_ssm)* np.ones(om_gw_ssm.shape)
#        
#        
#        param_estimate , _ = curve_fit(model,z,om_gw_ssm, bounds =(lower_bounds, upper_bounds),p0 = initial_guess, sigma = sigma )# curve_fit(func(xdata,p1,p2), xdata,ydata )
#        return om_gw_ssm, param_estimate
#    
#def fit_gw_ssm_2(z,params,model):
##    z,zp,Om_peak,zb
#    
#        om_gw_ssm = s.power_gw_scaled(z, params, npt = [5000,1000])
#        
#        # curve_fit set up 
#        
#        lower_bounds = [0,0,0]
#        upper_bounds = [np.inf,np.inf,5] # change back if doing rp !!!
#        
#        
#        initial_guess = [20,max(om_gw_ssm),1]
#        sigma = max(om_gw_ssm)* np.ones(om_gw_ssm.shape)
#        
#        
#        param_estimate , _ = curve_fit(model,z,om_gw_ssm, bounds =(lower_bounds, upper_bounds),p0 = initial_guess, sigma = sigma )# curve_fit(func(xdata,p1,p2), xdata,ydata )
#        return om_gw_ssm, param_estimate
#
#
#def fit_gw_ssm_3(z,params,model):
##   z,zb,Om_peak,bp
#    
#        om_gw_ssm = s.power_gw_scaled(z, params, npt = [5000,1000])
#        
#        # curve_fit set up 
#        
#        lower_bounds = [0,0,1]
#        upper_bounds = [5,np.inf,np.inf] # change back if doing rp !!!
#        
#        
#        initial_guess = [1,max(om_gw_ssm),20]
#        sigma = max(om_gw_ssm)* np.ones(om_gw_ssm.shape)
#        
#        
#        param_estimate , _ = curve_fit(model,z,om_gw_ssm, bounds =(lower_bounds, upper_bounds),p0 = initial_guess, sigma = sigma )# curve_fit(func(xdata,p1,p2), xdata,ydata )
#        return om_gw_ssm, param_estimate





#%%
#diff parameterisation models 

# Bounds for parameterisation #1
##z,zp,Om_peak,rp
#
#lower_b1 = [0,0,0]
#upper_b1 = [np.inf,np.inf,1]
#
#initial_guess1 = [1,0,0.5]
#
#
## Bounds for parameterisation #2
##z,zp,Om_peak,zb
#lower_b2 = [0,0,0]
#upper_b2 = [np.inf,np.inf,5]
#
#initial_guess2 = [1,0,0.5]
#
## Bounds for parameterisation #3
##z,zb,Om_peak,bp
#lower_b3 = [0,0,1]
#upper_b3 = [5,np.inf,np.inf]
#
#initial_guess3 = [0.5,0,2]



# Bounds for parameterisation #4 bp rescaled 
#z,Om_peak,rb , constant zb
#lower_2p = [0,0]
#upper_2p = [np.inf,1]
#
#initial_guess_2p = [0,0.5]




#%%

#plt.figure()
#for i, alpha in enumerate(alphas): 
#    plt.plot(vws,param_estimates_rescaled[i,:,0], label=  alpha )
#plt.xlabel(r'$v_w$')
#plt.ylabel(r'$z_b$')
#
#plt.legend(title= r'$\alpha$')
#plt.tight_layout()

#
#
#plt.figure()
#for i, alpha in enumerate(alphas): 
#    plt.plot(vws,1/param_estimates_rescaled[i,:,2], label=  alpha )
#    
#    
#plt.xlabel(r'$v_w$')
#plt.ylabel(r'$\frac{1}{r_b} = b+p$')
#plt.legend(title= r'$\alpha$')
#plt.tight_layout()
#

#plt.figure()
#for i, alpha in enumerate(alphas): 
#    plt.plot(vws,param_estimates[i,:,0]*param_estimates[i,:,2], label=  alpha )
#    
#    
#plt.xlabel(r'$v_w$')
#plt.ylabel(r'$z_p$')
#plt.legend(title= r'$\alpha$')
#plt.tight_layout()

#

#%%
    
#    
#vw_list = [vws[55],vws[70],vws[72],vws[74],vws[75],vws[76],vws[77],vws[78],vws[79],vws[80]]
#alpha_list =  0.06 * np.ones(len(vw_list))
##    
##plt.figure()
##for i,alpha in enumerate(alphas):
##    
#b.plot_fluid_shells(vw_list,alpha_list)

#%%


#
#
#cs = 1/np.sqrt(3)
#polyfitter(vws-cs,alphas,param_estimates_rescaled[:,:,3,2])



#%%

#    
##%% 
## what I'd done before ( 1/rb)
#
#om_tilde_2 =  np.zeros([n_alpha, n_vw])
#plt.figure()
#for i,alpha in enumerate(alphas):
#
#    for j,vw in enumerate(vws): 
##                                           omp,vw,alpha,1/rb
#        om_tilde_2[i,j] = g.om_to_om_tilde(param_estimates_rescaled[i,j,1],vw,alpha, 1/param_estimates_rescaled[i,j,2])
#    plt.semilogy(vws,om_tilde_2[i,:], label=  alpha )
#
#    
#plt.ylabel(r'$\tilde{\Omega}_{gw}$')
#plt.title('previous 1/rb')
#get_plot_setup()




##%%
##
#plt.figure()
#plt.loglog(z,om_gw_ssm[2,1,:]/P)
#plt.loglog(z,om_gw_ssm[2,10,:]/P)
#plt.loglog(z,om_gw_ssm[2,20,:]/P)
#plt.loglog(z,om_gw_ssm[2,30,:]/P)
#plt.loglog(z,om_gw_ssm[2,40,:]/P)
#plt.loglog(z,om_gw_ssm[2,50,:]/P)
#plt.loglog(z,om_gw_ssm[2,60,:]/P)
#plt.loglog(z,om_gw_ssm[2,70,:]/P)
#plt.loglog(z,om_gw_ssm[2,80,:]/P)
#plt.loglog(z,om_gw_ssm[2,90,:]/P)
#plt.loglog(z,om_gw_ssm[2,99,:]/P)