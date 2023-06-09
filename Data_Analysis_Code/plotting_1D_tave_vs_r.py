############################################################################################################################################  
#THIS SCRIPT PLOTS THE TIME AVERAGE QUNANITY VS R AT A FIXED RADIUS
############################################################################################################################################ 
# Arguments:
# 1. case (Wedge8, Wedge9B, Wedge10)  2. start iteration  3. end iteration 
# 4-7. list of reading/loading quantities  8. compute Mdot?
# 9. read/load time files? (set to False if we want to generate quantity files only)  10. succinct output?
############################################################################################################################################ 
#**Example command 1 (plotting Mdot vs r): python plotting_1D_tave_vs_r.py Wedge8 2000 3568 rhovr None None None True True False

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
                                                                                                                   
from parameters import *                                                                                           
from util import *                     

Print_title("========================== Starting plotting_1D_tave_vs_r.py =============================") 

########################################################################################################         
# input arguments                                                                                                
########################################################################################################       
Input_Arguments(argv)

case = argv[1]       
start = int(argv[2])                                                                                            
end = int(argv[3])                                                                                              
quant_list = [argv[4], argv[5], argv[6], argv[7]] # the last 3 are optional 
comp_Mdot = argv[8]
read_time = argv[9]
succinct = argv[10]

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
# Getting the theta array:                                                                            
# Any hist file can be used to generate these:                                                                   
########################################################################################################         
filename = "hist_"+str(start).zfill(5)+".npz"                                                                    
if succinct == "False": Print_subtitle("Getting Radius and Theta arrays, using file:", filename)      
data = np.load(dir+'/'+filename)                                                                                 

r  = Get_All_1D('radius', data, -1, dir, succinct)
th = Get_All_1D('theta', data, -1, dir, succinct)                                                                      
########### Get the differential theta array ###########################################################
dth = list(np.diff(th))                                                                                                         
dth.append(dth[-1])


########################################################################################################                         
# Getting the time array and quantity arrary:
########################################################################################################
quant_list = list(filter(lambda a: a != "None", quant_list))
radius = {'1': 120, '2':140, '3':160, '4':180}  

for idx, item in enumerate(quant_list):
    if comp_Mdot == 'True':
        quant_list[idx:idx] = ['Mdot']
        quant_list.remove('rhovr')

if succinct == "False": Print_subtitle("The datasets we want to read in are:", quant_list)

t = []                                                                                                                   
quant_data = []                                                                                                          
save_start_list = []                                                                                                                                                                                        
read_start_list = [] 
for item in quant_list:                                                                                                  
    quant_data.append([]) 
    save_start_list.append([])                                                                                                                                                                              
    read_start_list.append([])

file_exist = 0
for idx, quant in enumerate(quant_list):
    Print_subtitle("Compute " + str(quant) + "!!! First check saved files")     
    filenames = checkpoint_path + "T_ave_vs_r"+str(quant)+"_*" 
    t_filenames_pre = checkpoint_path + "Tave_vs_r_time"+str(quant)+"_"

    quant_data[idx], save_start_list[idx], read_start_list[idx] = Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct)


if read_time == "True": 
    t, time_save_start, time_read_start = Check_Load_Files(t_filenames_pre+"*", t_filenames_pre, file_exist, start, end, False, read_time, succinct)

start = min(read_start_list)    
if succinct == "False":                                                                                           
    Print_subtitle("Data structure after loading saved files:")                                                  
    Print_text("Time array:", np.shape(t)[0])                                                                    
    Print_text("The list of quantities and the shape of their data set:")                                                                                                                                   
    Print_text( [ (q, np.shape(quant_data[quant_list.index(q)])) for q in quant_list] )                                                                                                                     
    Print_text("Save start and read start lists are:", save_start_list, read_start_list)
    Print_text("(!) Starting iteration for new computation:", start)                                             
    Print_text("(!) Expected to read and plot data until iteration:", end)

save_step = 10
save_quant_file_pre = checkpoint_path + "T_ave_vs_r"+str(quant) + "_"
save_time_file_pre = checkpoint_path + "Tave_vs_r_time"+str(quant) + "_"

pbar = tqdm(total=end-start+1)
for iter in range(start, end):                                                                                   
    filename = "hist_"+str(iter).zfill(5)+".npz"                                                                 
    if succinct == "False":
        Print_subtitle("Reading file:", filename)
    else: 
        pbar.update(n=1)  
    data = np.load(dir+'/'+filename)
    idx_Mdot = -1
    #############################################
    #  Check the quantity data is 1D or 2D    
    #############################################
    for idx, quant in enumerate(quant_list):
        if quant != 'None' and iter >= read_start_list[quant_list.index(quant)]:
                if quant == 'Mdot':                    
                    idx_Mdot = idx
                    fac = 2 * np.pi * Mdot_to_cgs/Mdot_Edd_cgs
                    zip_data = list(zip(*data['rhovr']))
                    integrals = []
                    for i in range(len(zip_data)):
                       integral = np.sum([rhovr*np.sin(theta)*dtheta for rhovr, theta, dtheta in zip(zip_data[i], th, dth)])
                       integrals.append(integral*fac*r[i]**2)                   
                    #################### Mdot = 2pi * r^2 * Sum( rhovr[theta] * sintheta * dtheta ) ################
                    if np.shape( quant_data[idx] )[0] == 0: 
                        quant_data[idx].append(integrals)
                    else:
                        quant_data[idx] = np.vstack([quant_data[idx], integrals]) 
                else:
                    if np.shape( quant_data[idx] )[0] == 0: 
                        quant_data[idx].append(list(data[quant]))            
                    else:
                        quant_data[idx] = np.vstack([quant_data[idx], list(data[quant])])

    # Get time after Mdot data in case saved time files exists
    if read_time == "True" and iter >= time_read_start:                                                                                                                                                 
        if len(t) == 0:
            t.append(Get_Time(data, dir))  
        else:
            t = np.append(t, Get_Time(data, dir))
    
    if succinct == "False": Print_subsubtitle("Finish Reading!", "Time array length:", np.shape(t)[0], " Quantity arrays lengths",[[quant_list[i], np.shape(quant_data[i])] for i in range(len(quant_list))])  

########################################################################################################                               
# Save files since appending above is too slow!!!                                      
########################################################################################################
    if iter > start and (iter%save_step == 0 or iter == end-1):
        Save_Files(save_step, iter, save_start_list[idx], quant_data[idx_Mdot], save_quant_file_pre, succinct)
        
        if read_time == "True": 
            Save_Files(save_step, iter, time_save_start, t, save_time_file_pre, succinct)
                 
if succinct == "False": Print_subsubtitle(np.shape(quant_data), np.shape(quant_data[0]), quant_list)



########################################################################################################    
# Time average the quantities
########################################################################################################
ave_quant_data = []                                                                                             
for _ in range(len(quant_list)):                                                                            
    ave_quant_data.append([]) 

for idx, lists in enumerate(quant_data):
    if len(lists) > 0:
        ave_quant_data[idx] = list(np.mean(lists, axis = 0) )

ave_quant_data = [q for q in ave_quant_data if q != []]   


Print_title(sp30, "Start Plotting!!!!", sp30)  
########################################################################################################    
# Plotting and save
########################################################################################################
plot_dic = {"PB":['black', r'$P_B$', P_to_cgs, True],                                                                                                                                                               
            "pg":['blue', r'$P_{gas}$', P_to_cgs, True],                                                                                                                                                            
            "rhovt2":['red', r'$\rho v^2_{\theta}$', P_to_cgs, True],                                                                                                                                               
            "ST_PB1":['teal', r'$P_B^r$', P_to_cgs/2.0, True],                                                                                                                                                      
            "ST_PB2":['purple', r'$P_B^{\theta}$', P_to_cgs/2.0, True],                                                                                                                                             
            "ST_PB3":['steelblue', r'$P_B^{\phi}$', P_to_cgs/2.0, True],                                                                                                                                            
            "PB1":['teal', r'$P_B^r$', P_to_cgs, True],                                                                                                                                                             
            "PB2":['purple', r'$P_B^{\theta}$', P_to_cgs, True],                                                                                                                                                    
            "PB3":['steelblue', r'$P_B^{\phi}$', P_to_cgs, True],                                                                                                                                                   
            "Ek1":['peru', r'$E_k^r$', P_to_cgs, True],                                                                                                                                                             
            "Ek2":['darkorange', r'$E_k^{\theta}$', P_to_cgs, True],                                                                                                                                                
            "Ek3":['chocolate', r'$E_k^{\phi}$', P_to_cgs, True],                                                                                                                                                   
            "PB_total":['maroon', r'$P_B^{total}$', P_to_cgs, True],                                                                                                                                                
            "Er":['crimson', r'$P_{rad}$', Er_to_Pr*P_to_cgs, True],                                                                                                                                         
            "Mdot":['teal', r'$\hline{\dot{M}}/\dot{M}_{Edd}$', 1.0, False]
           } 
time_start = str(int(t[0]*time_to_sec))
time_end = str(int(t[-1]*time_to_sec))

fig, ax1 = plt.subplots(1, 1) 
ax1.set_title(str(case)+': Shell Averaged Accretion Rate From t = '+time_start+' to '+time_end, fontdict = font2)
xlabel = r'$r/r_g$'
ylabel = r'$\overline{\dot{M}}/\dot{M}_{Edd}$'
plt.xlim(100,2005)


for q_idx, item in enumerate(quant_list):
    item_pref = item.split("_")[0]
    color, label_q, q_fac, log_q = plot_dic[item_pref]  
    label_q = None
    Plotting_1D(ax1, r, 1.0, False, ave_quant_data[q_idx], q_fac, log_q, color, xlabel, ylabel, label_q, font)  
    if item == "Mdot":
       d0 = ave_quant_data[idx][0]
       format_d0 = "{:.2f}".format(d0) 
       ax1.text(r[0]*1.15, -5.5, r'$\overline{\dot{M}}/\dot{M}_{Edd}(r = r_{in})=$' + str(format_d0), rotation = 90)    



savename = save_path + 'Time_Averaged_'+str(quant_list[0])+str(case)+'vs_r'+time_start+'_'+time_end+'.png'
plt.savefig(savename,format='png',dpi=300)



