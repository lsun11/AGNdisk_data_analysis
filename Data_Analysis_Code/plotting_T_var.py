#THIS SCRIPT PLOTS THE TEMPERATURE VARIATIONS VS TIME AND RADIUS
#Location can be chosen at: midplane, theta = 80deg, theta = 70deg, and photosphere
#Example command: python plotting_T_var.py Wedge8 Er 3000 3050 photosphere rho sigma -1 1100 False

import numpy as np                                                                                                                                                                                                  
import os                                                                                                                                                                                                           
import sys                                                                                                                                                                                                          
import matplotlib.pyplot as plt                                                                                                                                                                                     
#import glob                                                                                                                                                                                                        
from sys import argv                                                                                                                                                                                                
myfonts = "Times New Roman"                                                                                                                                                                                         
plt.rcParams['font.family'] = "sans-serif"                                                                                                                                                                          
plt.rcParams['font.sans-serif'] = myfonts                                                                                                                                                                           
from tqdm import tqdm                                                                                                                                                                                               
from tqdm.gui import tqdm_gui

from parameters import *                                                                                                                                                                                            
from util import *     

Print_title("========================== Starting plotting_T_var.py =============================")


########################################################################################################                 
# input arguments                                                                                                        
########################################################################################################
Input_Arguments(argv)


case = argv[1]                                                                                                                                                                                                                    
quant = argv[2]   #Er for T
start = int(argv[3])
end = int(argv[4])                                                                                                                                                                                                                                  
deg = argv[5]    #midplane/80/70/photosphere
quant2 = argv[6] #rho for computing photosphere
quant3 = argv[7] #kappa for computing photosphere
ph_mode = int(argv[8])
rad_cut = int(argv[9])
succinct = argv[10]
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
dir = Data_dir + 'DATA' + str(case)                                                                                                                                                                                                            
######################################################################################################## 

########################################################################################################                 
#  set degree indices  (Currenly hard-coded since the size of array if fixed to be 2048)
########################################################################################################
if deg == "midplane":                                                                                                    
    deg_idx = 1023                                                                                                       
elif deg == "80":                                                                                                        
    deg_idx = 1023*8//9                                                                                                  
elif deg == "70":                                                                                                        
    deg_idx = 1023*7//9 
########################################################################################################                 

########################################################################################################  
# Getting the radius and theta array:
# Any hist file can be used to generate these:
########################################################################################################  
filename = "hist_"+str(start).zfill(5)+".npz"
Print_text(filename)                                                                                                      
data = np.load(dir+'/'+filename)                         

r = Get_All_1D('radius', data, rad_cut, dir, succinct)   
if deg == 'photosphere':
    th = Get_All_1D('theta', data, -1, dir, succinct)
########################################################################################################  


######################################################################################################## 
print("start reading time and data! data is", quant)
t          = []                                                                                                                                                                                       
quant_data = []
file_exist = 0
########################################################################################################                 
# If we plot over photosphere, we check if there are previously saved files                                                                       
# If so we load them and continue. 
######################################################################################################## 
if deg == "photosphere":  
# First check if we have saved files that includes the current iter 
    print("Plotting at photosphere!!! First check saved files")
    filenames = checkpoint_path + "T_var_Er*"
    t_filenames_pre = checkpoint_path + "T_var_time_"
    
    quant_data, save_start, q_start = Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, True, succinct)
    t, time_save_start, t_start = Check_Load_Files(t_filenames_pre+'*', t_filenames_pre, file_exist, start, end, False, True, succinct)

print(t)
start = np.min([q_start, t_start])
######################################################################################################## 




########################################################################################################
# Iterate from start to end to compute and save time and quant
######################################################################################################## 
if succinct == "True": pbar = tqdm(total=end-start, desc='Files read')
for iter in range(start, end):
        filename = "hist_"+str(iter).zfill(5)+".npz"                                                                                                                                                                              
        if succinct == "True": 
            pbar.update(n=1)
        else:
            Print_subtitle(filename)     
        data = np.load(dir+'/'+filename)      

        if iter >= t_start:
            if len(t) > 0:
                t = np.append(t, Get_Time(data, dir))
            else:
                t.append(Get_Time(data, dir))


        if deg != "photosphere":
            if succinct == "False": Print_text("We are plotting quantities at "+str(deg)+"!!!") 
            quant_data.append(data[quant][deg_idx][0:rad_cut]) #midplane

        else:   #integrate rho*kappa to get photosphere
            if succinct == "False": Print_text("We are plotting quantities at the photosphere!!!")

            r_data    = r[0:rad_cut]   
            kappa_zip = list(zip(*data[quant3]))
            kappa_data = kappa_zip[0:rad_cut]      
            th_data   = th        
                      
            # initialize an array for the indices of photosphere theta
            iphot = []
            hemi = 'upper'

            for it_r in range(len(r_data)):
                if ph_mode < 0:
                    iphot.append(Get_Photosphere(th_data, kappa_data[it_r], r_data[it_r], hemi))
                else:
                    iphot.append(Get_Photosphere2(th_data, kappa_data[it_r], r_data[it_r], mode, hemi)) 
 
            # This line below is SUPER SLOW! 
            if np.shape(quant_data)[0] == 0:
                phot_data = []
                if succinct == "False": pbar2 = tqdm(total=len(iphot), desc='Reading photosphere thetas in this file')
                for i in range(len(iphot)):                                                                                         
                    if succinct == "False": pbar2.update(n=1)                                                                                               
                    phot_data.append(data[quant][iphot[i]][i])    
                quant_data.append(phot_data)
                #quant_data.append([data[quant][iphot[i]][i] for i in range(len(iphot))]) 
                
                if succinct == "False": Print_subsubtitle("The structure of the read data:", np.shape(quant_data))

            else:
                if succinct == "False": pbar3 = tqdm(total=len(iphot), desc='Reading photosphere thetas in this file')
                phot_data = []
                for i in range(len(iphot)):
                    if succinct == "False": pbar3.update(n=1) 
                    phot_data.append(data[quant][iphot[i]][i])
                quant_data = np.vstack([quant_data,  phot_data])
                #quant_data = np.vstack([quant_data, [data[quant][iphot[i]][i] for i in range(len(iphot))]])
                #pbar2.update(n=1)

######################################################################################################## 
# Save files since appending above is too slow!!!
######################################################################################################## 
            save_step = 10                       
            if iter > start and (iter%save_step == 0 or iter == end-1):             
                save_quant_file_pre = checkpoint_path + "T_var_Er"
                Save_Files(save_step, iter, save_start, quant_data, save_quant_file_pre, succinct)

                save_time_file_pre = checkpoint_path + "T_var_time"
                Save_Files(save_step, iter, save_start,  t, save_time_file_pre, succinct) 

########################################################################################################
# Set file_exist to 0. This is necessary if the script is called multiple times at once 
########################################################################################################
file_exist = 0

######################################################################################################## 
# Converting Data
######################################################################################################## 
if deg == "photosphere":
    print("Converting quant data to arrays")
    quant_data = np.array(quant_data)


T_data = [q**0.25 for q in quant_data]

T_ave = [np.average(T) for T in T_data]
T_var = []

for i in range(len(T_data)):
    T_var.append([T - T_ave[i] for T in T_data[i]])


T_var = list(zip(*T_var))
######################################################################################################## 


######################################################################################################## 
# Plotting and save figure
######################################################################################################## 

x = t
y = r
q = T_var
x_fac = time_to_sec
y_fac = 1.0
q_fac = 1.0
log_x = False
log_y = True
log_q = False
color = 'RdBu_r'
xlim = None
ylim = None
v_max = None
v_min = None
title = r'$(\delta T_{rad}/T_{rad})/Scale$'
cbar_title = r'$(\delta T_{rad}/T_{rad})/Scale$'  
xlabel = 'Time (s)'
ylabel = '$\log_{10}(r/r_g)$'


time_start = int(t[0]*time_to_sec)
savename = save_path+'photosphere_rad_temperature_var_'+str(case)+'_'+str(deg)+'_t_'+str(time_start)+'.png'
font_dict = font

Plotting_Mesh(x,y,q,x_fac,y_fac,q_fac,log_x,log_y,log_q, color, xlim, ylim, v_max, v_min, title, xlabel, ylabel, cbar_title, font_dict)
plt.savefig(savename,format='png',dpi=300)

'''
Time, Radius = np.meshgrid(t, r)
Time *= time_to_sec                   #convert to second
time_start = int(t[0]*time_to_sec)
print(np.shape(Time), np.shape(Radius))
 

fig, ax1 = plt.subplots(1, 1)
im0 = ax1.pcolormesh(Time, np.log10(Radius), T_var, cmap='RdBu_r', shading='auto')
fig.colorbar(im0, ax=ax1)  
ax1.set_title(r'$(\delta T_{rad}/T_{rad})/Scale$')                                                                                                                                                                            
ax1.set_xlabel('Time (s)') 
ax1.set_ylabel('$\log_{10}(r/r_g)$')                                                                                                                                                                                          
#ax1.set_aspect(0.7)     
plt.savefig(save_path+'photosphere_rad_temperature_var_'+str(case)+'_'+str(deg)+'_t_'+str(time_start)+'.png'.format(i),format='png',dpi=300)
'''

######################################################################################################## 
print("========================== End plotting_T_var.py =============================") 
######################################################################################################## 

