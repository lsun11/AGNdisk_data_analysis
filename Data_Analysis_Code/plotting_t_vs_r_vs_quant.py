#THIS SCRIPT PLOTS TIME(X), RADIUS(Y), VS QUANTUTY, AT A GIVEN RADIUS
# Example1 command python plotting_t_vs_r_vs_quant.py Wedge8_2 surface_density None 90.0 3000 3500 False
# Example2 command python plotting_t_vs_r_vs_quant.py Wedge8_2 kappa rho 90.0 3000 3500 False

print("========================== Starting plotting_t_vs_r_vs_quant.py =============================")
import numpy as np
import os
import matplotlib.pyplot as plt
import sys
from sys import argv
myfonts = "Times New Roman"                                                                                                                                                                                                            
plt.rcParams['font.family'] = "sans-serif"                                                                                                                                                                                             
plt.rcParams['font.sans-serif'] = myfonts                                                                                                                                                                                              
                                                                                                                                                                                                                                       
from parameters import *                                                                                            
from util import *
                                                                                                                                                                                                                                       
font = {'family': 'sans-serif',                                                                                                                                                                                                        
        'color':  'black',                                                                                                                                                                                                             
        'weight': 'normal',                                                                                                                                                                                                            
        'size': 16,                                                                                                                                                                                                                    
        }                                                                                                                                                                                                                              
    
########################################################################################################            
# input arguments                                                                                                   
########################################################################################################            
for arg in sys.argv[:]:                                                                                             
   if arg.endswith('.txt'):                                                                                         
       with open(arg[:]) as f:                                                                                      
           file_args = f.read().splitlines()                                                                        
           argv[:] = [arg for arg in file_args if arg[0] != "#"]                                                    
   else:                                                                                                            
       argv.append(arg)

case = argv[1]
quant1 = argv[2]
quant2 = argv[3]        # optional
theta_select = argv[4]  # in degree, from 0 to 180
start = int(argv[5])
end   = int(argv[6])
succinct = argv[7]
######################################################################################################## 


########################################################################################################            
#  set paths based on cases                                                                                         
########################################################################################################
if str(case).find('Wedge8') == 0:                                                                                     
    print('case: Wedge8')                                                                                         
    checkpoint_path = Checkpoint_dir  + 'Wedge8/'                                                                     
    save_path = Plot_dir + 'Wedge8/'                                                                              
elif str(case).find('Wedge9B') == 0:                                                                                  
    print('case: Wedge9B')                                                                                        
    checkpoint_path = Checkpoint_dir  + 'Wedge9B/'                                                                    
    save_path = Plot_dir + 'Wedge9B/'                                                                             
elif str(case).find('Wedge9Res') == 0:                                                                                
    print('case: Wedge9Res')                                                                                      
    checkpoint_path = Checkpoint_dir  + 'Wedge9Res/'                                                                  
    save_path = Plot_dir + 'Wedge9Res/'                                                                           
elif str(case).find('Wedge10') == 0:                                                                                  
    print('case: Wedge10')                                                                                        
    checkpoint_path = Checkpoint_dir  + 'Wedge10/'                                                                    
    save_path = Plot_dir + 'Wedge10/' 
dir = Data_dir + 'DATA'+str(case)
########################################################################################################


########################################################################################################            
# Getting the radius and theta array:                                                                               
# Any hist file can be used to generate these:                                                                      
########################################################################################################            
filename = "hist_"+str(start).zfill(5)+".npz"                                                                       
print(filename)                                                                                                     
data = np.load(dir+'/'+filename)                                                                                    
                                                                                                                    
r = Get_All_1D('radius', data, -1, dir, succinct) 
th = Get_All_1D('theta', data, -1, dir, succinct)                                                                 

########################################################################################################                         
# Find the index of theta_select
########################################################################################################
shift_theta = [np.abs(theta*rad_to_deg - float(theta_select)) for theta in th]   
theta_idx = shift_theta.index(min(shift_theta))   
########################################################################################################

t 	   = []
if str(case)[-2:] == '_2':
    print('case: _2')
    idx = {"surface_density":6, 'rho':7, 'Er':8, 'kappa':9, 'Fr1':10, 'Fr2':11, 'B1_r1':12, 'B1_r2':13, 'B1_r3':14, 'B1_r4':15, 'B2_r1':16, 'B2_r2':17, 'B2_r3':18, 'B2_r4':19, 'B3_r1':20, 'B3_r2':21, 'B3_r3':22, 'B3_r4':23, 'PB_r1':24, 'PB_r2':25, 'PB_r3':26, 'PB_r4':27, 'pg_r1':28, 'pg_r2':29, 'pg_r3':30, 'pg_r4':31, 'rhovr':32, 'lambda_r1':33, 'lambda_r2':34, 'lambda_r3':35, 'lambda_r4':36}
else:
    print('case: no _2')
    idx = {"surface_density":6, 'rho1':7, 'rho2':8, 'rho3':9, 'rho4':10, 'Er1':11, 'Er2':12, 'Er3':13, 'Er4':14, 'kappa1':15, 'kappa2':16, 'kappa3':17, 'kappa4':18}

quant1_data = []
quant2_data = []

for iter in range(start, end):
    filename = "hist_"+str(iter).zfill(5)+".npz"
    print(filename)
    data = np.load(dir+'/'+filename)  
    t.append(Get_Time(data, dir)) 

    if len(np.shape(data[quant1])) == 1:     # 1D data (integrated for all theta)
        quant1_data.append( data[quant1] )
    else:                                  # 2D data
        quant1_data.append( data[quant1][theta_idx] )  
    if quant2 != "None":
        if len(np.shape(data[quant1])) == 1: # 1D data (integrated for all theta)  
            quant2_data.append( data[quant2] )
        else:                              # 2D data                      
            quant2_data.append( data[quant2][theta_idx] ) 

print("Done reading data!")
print(np.shape(quant1_data), np.shape(quant2_data))


if quant1 == "kappa" and quant2 == "rho":
    print("Computing opacity!")
    quant_data = [[kk/rr for kk, rr in zip(k, rho)] for k, rho in zip(quant1_data, quant2_data)]
else:
   quant_data = quant1_data


if str(str(quant1)[0:3]) == "sur":   #surface_density
    color = 'RdBu_r'
    pre_str = '\\log_{10}(\\Sigma) (g/cm^2)'
    cbar_str = r'$\log_{10}(\Sigma) (g/cm^2)$'
    fac = dens_to_cgs
    logscale = True 
elif str(str(quant1)[0:3]) == "rho":  
    color = 'BrBG_r'                                                                                             
    pre_str = '\\rho/\\rho_0'
    cbar_str = r'$\rho/\rho_0$'                                                                                   
    fac  = 1.0/dens_to_cgs                                                                                             
    logscale = True                                                                                              
elif str(str(quant1)[0:2]) == "Er":   
    color = 'PRGn_r'                                                                                             
    pre_str = 'E_r'
    cbar_str = r'$E_r$'                                                                                           
    fac = 2.77e5                                                                                                 
    logscale = True                                                                                              
elif str(str(quant1)[0:5]) == "kappa": 
    color = 'RdGy_r'  
    pre_str = '\\log_{10}(\\kappa) (cm^2/g)'                                                                                          
    cbar_str = r'$\log_{10}(\kappa)$' 
    fac = kappa_to_cgs
    logscale = True

title = r'${}$'.format(pre_str) 
print(title)
xlabel = r'$t (s)$'
ylabel = r'$\log_{10}(r/r_g)$'
y_lim=(np.log10(100), np.log10(500))

print(np.shape(quant1_data))
################## Logscale of q has been taken care of, set log_q to False here ###################################################################
ax1 = Plotting_Mesh_YX(t, r, quant1_data, time_to_sec, 1.0, fac, False, True, logscale, color, None, y_lim, None, None, title, xlabel, ylabel, cbar_str, font)

save_str = '_' + str(theta_select) + '_'  if len(np.shape(data[quant1])) == 2 else '_'
print (save_str)

savename = save_path + str(quant1) + '_vs_t_and_r_' + str(case) + str(save_str) + str(start) + '_' +str(end) + '_' + 'plot.png' 
plt.savefig(savename,format='png',dpi=300)
#####################################################################################################################################################

'''
fig, ax1 = plt.subplots(1, 1)
print(np.shape(Time), np.shape(Theta), np.shape(quant_data))
im0 = ax1.pcolormesh(Time, Theta, quant_data, cmap=color, shading='auto')
fig.colorbar(im0, ax=ax1) 

phot_th_upper = [th[it]*rad_to_deg for it in iphot_upper_total]
phot_th_lower = [th[it]*rad_to_deg for it in iphot_lower_total]

ax1.plot(t_phys, phot_th_upper, color = 'y')
ax1.plot(t_phys, phot_th_lower, color = 'y')

radius = {'1': 120, '2':140, '3':160, '4':180}

tail = str(quant1)[-1]

if tail in radius.keys():
    string_radius = str(radius[tail])
else:
    string_radius = str(int(radius_select))



ax1.set_title(r'${}$ at $r =$'.format(pre_str) + str(theta_select) +'$^o$' + ' from ' + str(start_time) + ' to ' + str(end_time), fontdict = font)
#ax1.set_ylabel('$\log_{10}(r/r_g)$')                                                                                                    
ax1.set_ylabel(r'$\theta(^o)$', fontdict = font)
ax1.set_xlabel(r'$t (s)$', fontdict = font)              
#ax1.set_aspect(0.003)                                                                                                                    
#plt.savefig( save_path + str(quant1) + 'vs_t_and_theta_' + str(case) + '_' + string_radius + '_' + str(start) + '_' +str(end) +'.png',format='png',dpi=300)
'''

print("========================== Done plotting_t_vs_th_vs_quant.py =============================")
