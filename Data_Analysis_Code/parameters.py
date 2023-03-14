# This script has all the constants, unit factors, and common directories used by almost every program.  

c = 29979245800.0

rad_to_deg = 57.2958
time_to_sec = 495.0766 
r_to_cgs = 1.4842e13
v_to_cgs = c/8.05e3
dens_to_cgs = 5.0e-10
kappa_to_cgs = 1.0/(1.4842e13*5.0e-10)
P_to_cgs = dens_to_cgs * v_to_cgs**2
Mdot_to_cgs = dens_to_cgs * r_to_cgs**2 * v_to_cgs  
#kappa_to_cgs = dens_to_cgs * r_to_cgs
# The kappa vaue below is probably not right. I guess it using AGNIron's data
# 5.0655e5*(1.4842e13*5e-10)/(1e-8*2*7.42e13), where 5.0655e5 is AGNIron's kappa_to_cgs
kappaes_code = 2523.14    # electron scattering opacity in code units (~0.34 cm^2/g * kappa_to_cgs )
Er_to_Pr = 109.194/3

M_BH = 1.0e8
Mdot_Edd_cgs = 1.4e18*M_BH                         # in g/s

# Dictories to save plots and checkpoints 
Data_dir             = '/home/lunan/mnt/AGNWedge8_product/'
Plot_dir             = Data_dir + 'PLOTS/'
Checkpoint_dir       = Data_dir + 'Checkpoint/'

# Standard font
font = {'family': 'sans-serif',                                                                                                                                                                                     
        'color':  'black',                                                                                                                                                                                          
        'weight': 'normal',                                                                                                                                                                                         
        'size': 16,                                                                                                                                                                                                 
        } 

# Smaller font size
font2 = {'family': 'sans-serif',                                                                                       
        'color':  'black',                                                                                            
        'weight': 'normal',                                                                                           
        'size': 14,                                                                                                   
        }