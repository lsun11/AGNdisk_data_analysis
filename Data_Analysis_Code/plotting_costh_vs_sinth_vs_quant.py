############################################################################################################################################ 
#THIS SCRIPT PLOTS R*SIN(THETA)(X), R*COS(THETA)(Y), VS QUANTITY, AT A GIVEN TIME
############################################################################################################################################ 
# Example command python plotting_costh_vs_sinth_vs_quant.py Wedge8_2 rho sigma_p 2000 Y


import numpy as np
import os
import matplotlib.pyplot as plt
import sys
from sys import argv
myfonts = "Times New Roman"                                                                                                                                                                                                            
plt.rcParams['font.family'] = "sans-serif"                                                                                                                                                                                             
plt.rcParams['font.sans-serif'] = myfonts                                                                                                                                                                                              
from tqdm import tqdm

                                                                                                                                                                                                                                       
from parameters import *                                                                                            
from util import *
                                                                                                                                                                                                                                       
Print_title("========================== Starting plotting_costh_vs_sinth_vs_quant.py =============================") 
    
########################################################################################################            
# input arguments                                                                                                   
########################################################################################################            
Input_Arguments(argv)

case = argv[1]
quant1 = argv[2]       # quantity to plot
quant2 = argv[3]       # aux quantitym can be None
iter = argv[4]   # select iteration(time) to make the plot
plot_phot = argv[5]     # "Y" or "N"
#start = int(argv[6])
#end   = int(argv[7])
######################################################################################################## 


########################################################################################################            
#  set paths based on cases                                                                                         
########################################################################################################
if str(case).find('Wedge8') == 0: 
    Print_subtitle('We are analyzing the case: Wedge8')                                                           
    checkpoint_path = Checkpoint_dir  + 'Wedge8/'                                                                     
    save_path = Plot_dir + 'Wedge8/' 
elif str(case).find('Wedge9B') == 0: 
    Print_title('We are analyzing the case: Wedge9B')                                                             
    checkpoint_path = Checkpoint_dir  + 'Wedge9B/'                                                                    
    save_path = Plot_dir + 'Wedge9B/'
elif str(case).find('Wedge9Res') == 0:                                                                                
    Print_title('We are analyzing the case: Wedge9Res')                                                           
    checkpoint_path = Checkpoint_dir  + 'Wedge9Res/'                                                                  
    save_path = Plot_dir + 'Wedge9Res/'           
elif str(case).find('Wedge10') == 0:                                                                                  
    Print_title('We are analyzing the case: Wedge10')                      
    checkpoint_path = Checkpoint_dir  + 'Wedge10/'                                                                    
    save_path = Plot_dir + 'Wedge10/' 
dir = Data_dir + 'DATA' + str(case)


########################################################################################################            
# Getting the radius and theta array:                                                                               
# Any hist file can be used to generate these:                                                                      
########################################################################################################            
filename = "hist_"+str(iter).zfill(5)+".npz"                                                                       
Print_subtitle("Getting Radius and Theta arrays, using file:", filename) 
data = np.load(dir+'/'+filename)                                                                                    
                                                                                                                    
r = Get_All_1D('radius', data, -1, dir)                   
th = Get_All_1D('theta', data, -1, dir)                                                                 

        
########################################################################################################
# Getting the time and quantities arrays:                                                                         
########################################################################################################
t 	    = []
quant1_data = []
quant2_data = []


filename = "hist_"+str(iter).zfill(5)+".npz"
Print_subtitle("Getting Time array using file:", filename) 
t.append(Get_Time(data, dir)) 
t[0] *= time_to_sec

data = np.load(dir+'/'+filename)

quant1_data = data[quant1]
if quant2 != "None":
    quant2_data = data[quant2]

iphot_upper_total = []                                                                                                                                                                                              
iphot_lower_total = []          
if plot_phot == "Y":     # make sure your data[quant2] is kappa (or kappa*rho in reality) here!!!!
    Print_subsubtitle("Compute Photosphere!")
    pbar = tqdm(total=len(r))
    for radius_select in r:
        pbar.update(n=1)
        
        shift_rad = [np.abs(rad - radius_select) for rad in r]                                                                                                                                                 
        radius_idx = shift_rad.index(min(shift_rad))         
        quant1_ph_data = [ q[radius_idx] for q in data[quant1]]
        quant2_ph_data = [ q[radius_idx] for q in data[quant2]]
   
        #iphot_upper = Get_Photosphere(th, quant2_ph_data, float(radius_select), 'upper')                                                                                                                          
        #iphot_lower = Get_Photosphere(th, quant2_ph_data, float(radius_select), 'lower')

        iphot_upper = Get_Photosphere2(th, quant2_ph_data, quant1_ph_data, float(radius_select), 'upper')
        iphot_lower = Get_Photosphere2(th, quant2_ph_data, quant1_ph_data, float(radius_select), 'lower')
        iphot_upper_total.append(iphot_upper)
        iphot_lower_total.append(iphot_lower)
        #print(iphot_upper_total)

Print_subtitle("Done reading data!")
Print_text("The shape of the datasets are:", np.shape(quant1_data), np.shape(quant2_data))


if str(str(quant1)[0:2]) == "PB":                                                                                
    color = 'RdGy_r'                                                                                             
    pre_str = 'P_B'                                                                                              
    cbar_str = r'$P_B$'                                                                                          
    fac = 2.77e5                                                                                                 
    logscale = True                                                                                              
elif str(str(quant1)[0:2]) == "pg":                                                                              
    color = 'PuOr_r'                                                                                             
    pre_str = 'P_gas'                                                                                            
    cbar_str = r'$P_{gas}$'                                                                                      
    fac = 2.77e5                                                                                                 
    logscale = True                                                                                              
elif str(str(quant1)[0:3]) == "rho":  # no _2                                                                    
    color = 'BrBG_r'                                                                                             
    pre_str = 'log_{10}(\\rho/\\rho_0)'                                                                          
    cbar_str = r'$log_{10}(\rho/\rho_0)$'                                                                        
    fac  = 1.0/dens_to_cgs                                                                                       
    logscale = True                                                                                              
elif str(str(quant1)[0:2]) == "Er":   # no _2                                                                    
    color = 'PRGn_r'                                                                                             
    pre_str = 'E_r'                                                                                              
    cbar_str = r'$E_r$'                                                                                          
    fac = 2.77e5                                                                                                 
    logscale = True                                                                                              
elif str(str(quant1)[0:2]) == "B1":                                                                              
    color = 'RdYlBu_r'                                                                                               
    pre_str = 'E_r'         
    cbar_str = r'$B_r$'                                                                                              
    fac = 2.64e3                                                                                                 
    logscale = False                                                                                             
elif str(str(quant1)[0:2]) == "B2":                                                                              
    color = 'RdYlGn_r'                                                                                               
    pre_str = 'B_{\\theta}'                                                                                      
    cbar_str = r'$B_{\theta}$'                                                      
    fac = 2.64e3                                                                                                 
    logscale = False                                                                                             
elif str(str(quant1)[0:2]) == "B3":                                                                              
    color = 'RdBu_r'                                                                                             
    pre_str = 'B_{\\phi}'                                                                                        
    cbar_str = r'$B_{\phi}$'                                                                                     
    fac = 2.64e3                                                                                                 
    logscale = False                                                                                             
elif str(str(quant1)[0:5]) == "kappa": # no _2                                                                   
    color = 'seismic'                                                                                            
    pre_str = '\\kappa'                                                                                          
    cbar_str = r'$\kappa$'                                                                                       
    fac = 8.05e3/5/3                                                                                             
    logscale = False 



#rcosth = [ [rad*np.cos(theta) for theta in th] for rad in r]
rcosth = [ [rad*np.cos(theta) for rad in r] for theta in th] 
#rsinth = [ [rad*np.sin(theta) for theta in th] for rad in r]
rsinth = [ [rad*np.sin(theta) for rad in r] for theta in th] 


radius = {'1': 120, '2':140, '3':160, '4':180}                                                                                                                                                                                   

title = 't =' + str(int(t[0]))
xlabel = r'$r\sin{\theta}/r_g$'
ylabel = r'$r\cos{\theta}/r_g$'
xlim = [100, 500]
ylim = [-500, 500]

################## Logscale of q has been taken care of, set log_q to False here ###################################################################
ax1 = Plotting_Mesh_MeshXY(rsinth, rcosth, quant1_data, 1.0, 1.0, fac, False, False, True, color, xlim, ylim, None, None, title, xlabel, ylabel, cbar_str, font)

savename = save_path + str(quant1) + 'vs_cost_and_sint_' + str(case) + '_' + str(int(t[0])) + '_' + 'plot.png' 
################## Plot the photosphere curve ######################################################################################################
if plot_phot == "Y":
    rcos_phot_th_upper = [ rad*np.cos(th[it]) for rad, it in zip(r,iphot_upper_total)]                                                                                                                                                                  
    rsin_phot_th_upper = [ rad*np.sin(th[it]) for rad, it in zip(r,iphot_upper_total)]
    rcos_phot_th_lower = [ rad*np.cos(th[it]) for rad, it in zip(r,iphot_lower_total)]                                                                                                                              
    rsin_phot_th_lower = [ rad*np.sin(th[it]) for rad, it in zip(r,iphot_lower_total)] 
    
    Plotting_1D(ax1, rsin_phot_th_upper, 1.0, False, rcos_phot_th_upper, 1.0, False, 'y', None, None, None, font)
    Plotting_1D(ax1, rsin_phot_th_lower, 1.0, False, rcos_phot_th_lower, 1.0, False, 'y', None, None, None, font) 

    ####################### Add an extra part of the filename ###########################
    split_name = savename.split('_')
    split_name.insert(-1, 'photosphere')
    savename = ('_').join(split_name)
################## Save the plot ###################################################################################################################
plt.savefig(savename,format='png',dpi=300)
#####################################################################################################################################################


Print_title("========================== Done plotting_t_vs_th_vs_quant.py =============================")
