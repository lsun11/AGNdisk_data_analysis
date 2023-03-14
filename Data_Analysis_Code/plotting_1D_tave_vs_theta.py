############################################################################################################################################ 
#THIS SCRIPT PLOTS THE TIME AVERAGE QUNANITY VS THETA (vs_theta) OR VS TIME (vs_t) AT A FIXED RADIUS
# Careful!!!!! Quantiy labelled with _r1 to _r4 corresponds to radii 120 140 160 180. Make sure they match with the radius you select!!!!!!
############################################################################################################################################
# Arguments:
# 1. case (Wedge8, Wedge9B, Wedge10)  2. start iteration  3. end iteration 
# 4. select radius in rg  5-8. list of reading/loading quantities  9. plotting against which quantity (theta or time)   
# 10. read/load time files? (set to False if we want to generate quantity files only)  11. Verbose output? 
############################################################################################################################################ 
#Example command 1 (plotting different pressure at r = 120rg): python plotting_1D_tave_vs_theta.py Wedge8 2000 3658 120 PB_r1 pg_r1 rhovt2 None vs_theta True False
#Example command 2 (plotting different pressure at r = 180rg): python plotting_1D_tave_vs_theta.py Wedge8 2000 3658 180 PB_r4 pg_r4 rhovt2 None vs_theta True False
#Example command 3 (plotting diff B-field components): python plotting_1D_tave_vs_theta.py Wedge8_2 2000 3658 140 B1_r2 B2_r2 B3_r2 None vs_theta True False
#Example command 4 (plotting TOTAL Magnetic Pressure): python plotting_1D_tave_vs_theta.py Wedge8 2000 3658 120 PB_total pg_r1 Ek2 Er vs_theta True False

#THIS CAN ALSO PLOT THETA AVERAGE QUANTITY VS T AT A FIX RADIUS, BY CHANGING PLOTTING MODE TO vs_t
# Notice that vs_t code can use the same saved checkpoint files as vs_theta, no need to read everything again!
#Example command 1 (plotting different pressure): python plotting_1D_tave_vs_theta.py Wedge8 2000 3658 120 PB_r1 pg_r1 rhovt rho vs_t True False
#Example command 2 (plotting diff B-field components): python plotting_1D_tave_vs_theta.py Wedge8_2 2000 3658 140 B1_r2 B2_r2 B3_r2 None vs_t True False

import numpy as np                                                                                                 
import os                                                                                                          
import sys                                                                                                         
import matplotlib.pyplot as plt                                                                                    
#import glob                                                                                                       
from sys import argv                                                                                               
myfonts = "Times New Roman"                                                                                        
plt.rcParams['font.family'] = "sans-serif"                                                                         
plt.rcParams['font.sans-serif'] = myfonts                                                                          
                                                                                                                   
from parameters import *                                                                                           
from util import *                     

Print_title("========================== Starting plotting_1D_tave.py =============================")


########################################################################################################         
# input arguments                                                                                                
########################################################################################################       
Input_Arguments(argv)

case  = argv[1]       
start = int(argv[2])                                                                                            
end   = int(argv[3])                                                                                              
radius_select = argv[4]
quant_list = [argv[5], argv[6], argv[7], argv[8]] # the last 3 are optional 
x_axis     = argv[9]
read_time  = argv[10]
verbose    = argv[11]

if x_axis == "vs_theta":
    plot_string = "Theta"
    plot_axis = 0
elif x_axis == "vs_t":       
    plot_string = "Radius"   
    plot_axis = 1
else:                                                                                                                 
    Print_subtitle("Wrong Plotting Axis! Exit!")
    sys.exit()

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
if verbose == "False": Print_subtitle("Getting Radius and Theta arrays, using file:", filename)                                                                                                  
data = np.load(dir+'/'+filename)                                                                                 

r  = Get_All_1D('radius', data, -1, dir)
th = Get_All_1D('theta', data, -1, dir)                                                                      


########################################################################################################
# Check and Load files
########################################################################################################
# First convert quant_list to what we want
# 1. Remove all "None"s
# 2. If PB_total is included, add ST_PB1, ST_PB2, ST_PB3 to the list
######################################################################################################## 
quant_list = list(filter(lambda a: a != "None", quant_list))

PB_total_flag = 0
rhovt2_flag = 0 

for idx, item in enumerate(quant_list):
    if item == "PB_total":
        PB_total_flag = 1
        if '_2' in str(case):
            quant_list[idx:idx] = ['PB1', 'PB2', 'PB3']
        else:
            quant_list[idx:idx] = ['ST_PB1', 'ST_PB2', 'ST_PB3']
        quant_list.remove('PB_total')
    if item == "rho_vt2":                                                                                             
        rhovt2_flag = 1      
        quant_list[idx:idx] = ['rhovt', 'rho']                   
        quant_list.remove('rho_vt2')

if verbose == "False": Print_subtitle("The datasets we want to read in are:", quant_list)


t = []                                                                                                                   
quant_data = [] 
for item in quant_list:                                                                                         
    quant_data.append([]) 

file_exist = 0
for idx, quant in enumerate(quant_list):
    Print_subtitle("Compute " + str(quant) + "!!! First check saved files")
    filenames = checkpoint_path + "T_ave_vs_th_"+str(quant)+"_"+str(radius_select)+"_*"
    t_filenames_pre = checkpoint_path + "Tave_vs_theta_time_"+str(radius_select)+"__" 

    if idx != len(quant_list)-1 :
        quant_data[idx] = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, verbose))
    else:
        t, quant_data[idx], file_exist, save_start, start = Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, True, read_time, verbose)
        if verbose == "False":
            Print_subtitle("Data structure after loading saved files:")
            Print_text("Time array:", np.shape(t)[0])
            Print_text("Quantity and its data set:",[[quant_list[idx], np.shape(quant_data[idx])] for idx in range(len(quant_list))])
            Print_text("Starting iteration existed in the saved data:", save_start)
            Print_text("(!) Starting iteration for new computation:", start)
            Print_text("(!) Expected to read and plot data until iteration:", end)
save_step = 20
save_time_file_pre = checkpoint_path + "Tave_vs_theta_time_"+str(radius_select)+"_"

Print_title(sp20*5)
########################################################################################################                         
# Computing the new time array and quantity arrays:
########################################################################################################
Print_subtitle("Start Reading New Data to Arrays!!!")
radius = {'1': 120, '2':140, '3':160, '4':180}  
for iter in range(start, end):                                                                                   
    filename = "hist_"+str(iter).zfill(5)+".npz"
    if os.path.exists(dir+'/'+filename):
        Print_subsubtitle("Reading file:", filename)
        data = np.load(dir+'/'+filename)
        #######################################
        #  Check the quantity data is 1D or 2D    
        #######################################
        for idx, quant in enumerate(quant_list):
            if quant != 'None':
                if verbose == "False": Print_text("Reading quanity "+str(quant))
                if len(np.shape(data[quant])) ==  1:   # 1D
                    if file_exist == 0:
                        quant_data[idx].append(list(data[quant]))
                    else:
                        quant_data[idx] = np.vstack([quant_data[idx], list(data[quant])])
                else:                                 #  2D
                    shift_rad = [np.abs(rad - float(radius_select)) for rad in r]  
                    radius_idx = shift_rad.index(min(shift_rad))

                    if file_exist == 0:
                        quant_data[idx].append(list([q[radius_idx] for q in data[quant]])) 
                    else:
                        quant_data[idx] = np.vstack([quant_data[idx], list([q[radius_idx] for q in data[quant]])])

                    #string_radius = str(int(radius_select))

        # Get time after Mdot data in case saved time files exists
        if file_exist == 1:  
            t = np.append(t, Get_Time(data, dir))    
        else: 
            t.append(Get_Time(data, dir)) 
        if verbose == "False": Print_subsubtitle("Finish Reading!", "Time array length:", np.shape(t)[0], " Quantity arrays lengths",[[quant_list[i], np.shape(quant_data[i])] for i in range(len(quant_list))])

########################################################################################################       
# Save files since appending above is too slow!!!  
########################################################################################################      
        if iter > start and (iter%save_step == 0 or iter == end-1): 
            for idx, quant in enumerate(quant_list):
                Print_subtitle("Saving files At Iteration", iter, "quantity: ", quant)
                Save_Files(save_step, iter, save_start, start, end, quant_data[idx], checkpoint_path + "T_ave_vs_th_"+str(quant)+"_"+str(radius_select)+"_")
 
            if read_time == "True":
                Print_subtitle("Saving files At Iteration", iter, "quantity: Time")
                Save_Files(save_step, iter, save_start, start, end, t, save_time_file_pre) 


########################################################################################################
# for rhovt^2 quantity, need (rhovt)^2/rho
# Average the quantity and append to ave_quant_data first
# If we append non-time averaged array first, we encounter ragged array problem 
######################################################################################################## 
if rhovt2_flag == 1: 
    Print_subtitle("Compute rhovt^2!!!")
    data_rhovt = quant_data[quant_list.index('rhovt')]
    data_rho = quant_data[quant_list.index('rho')] 

    #print(np.shape(quant_data), np.shape(data_rhovt), np.shape(data_rho))   
    data_rhovt2 = [ [rhovt**2/rho for rhovt, rho in zip (list_rhovt, list_rho)] for list_rhovt, list_rho in zip(data_rhovt, data_rho)]

    #print([ (i, np.shape(quant_data[i])) for i in range(len(quant_data))])
    ave_data_rhovt2 = np.mean(data_rhovt2, axis = plot_axis)
    ########### After getting rhovt^2,  we want to skip rhovt and rho for plotting ################
    quant_data = np.delete(quant_data, [quant_list.index('rhovt'), quant_list.index('rho')], 0)
    #print([ (i, np.shape(quant_data[i])) for i in range(len(quant_data))])      

    #quant_data = np.vstack([quant_data, [data_rhovt2]])
    quant_list.remove('rhovt')                                                                              
    quant_list.remove('rho')                                                                                    
    quant_list.insert(len(quant_list),'rhovt2') 
    if verbose == "False": Print_subsubtitle("After computing rhovt^2, Quantity arrays lengths",[[quant_list[i], np.shape(quant_data[i])] for i in range(len(quant_list))])
 
########################################################################################################                                                                                                            
# Compute PB_total = sqrt(ST_PB1^2 + ST_PB2^2 + ST_PB3^2)/2 = sqrt(PB1^2 + PB2^2 + PB3^2)
# Same. Average the quantity and append to ave_quant_data first  
######################################################################################################## 
if PB_total_flag == 1:
    print("Compute PB_total!!!")
    if "_2" in str(case):
       data_PB1 = quant_data[quant_list.index('PB1')]                                                                            
       data_PB2 = quant_data[quant_list.index('PB2')]                                                                                   
       data_PB3 = quant_data[quant_list.index('PB3')]
       quant_data = np.delete(quant_data, [quant_list.index('PB1'), quant_list.index('PB2'),  quant_list.index('PB3')], 0)
       quant_list.remove('PB1')                                                                                              
       quant_list.remove('PB2')       
       quant_list.remove('PB3')
    else:
       data_PB1 = quant_data[quant_list.index('ST_PB1')] 
       data_PB2 = quant_data[quant_list.index('ST_PB2')]
       data_PB3 = quant_data[quant_list.index('ST_PB3')]
       quant_data = np.delete(quant_data, [quant_list.index('ST_PB1'), quant_list.index('ST_PB2'),  quant_list.index('ST_PB3')], 0)  
       quant_list.remove('ST_PB1')
       quant_list.remove('ST_PB2')
       quant_list.remove('ST_PB3')
    
    quant_list.insert(len(quant_list),'PB_total')

    data_PB_total = [ [np.sqrt(PB1**2 + PB2**2 + PB3**2) for (PB1, PB2, PB3) in zip(list1, list2, list3)] 
                    for (list1, list2, list3) in zip (data_PB1, data_PB2, data_PB3) ]    
    ave_data_PB_total = np.mean(data_PB_total, axis = plot_axis) 


########################################################################################################    
# Time average the quantities
########################################################################################################
ave_quant_data = []                                                                                             
for _ in range(len(quant_list)):                                                                            
    ave_quant_data.append([]) 

for idx, lists in enumerate(quant_data):
    if len(lists) > 0:
        ave_quant_data[idx] = list(np.mean(lists, axis = plot_axis) )
    if rhovt2_flag == 1:                                                                                              
        ave_quant_data[quant_list.index('rhovt2')] = ave_data_rhovt2                                                  
    if PB_total_flag == 1:                                                                                            
        ave_quant_data[quant_list.index('PB_total')] = ave_data_PB_total 

if verbose == "False":
    Print_subtitle("The final qunaities for plotting are:", quant_list)
    Print_subsubtitle("Checking the plotting quantity data structure:")
    Print_text("Time, theta, radius:", np.shape(t), np.shape(th), np.shape(r))
    Print_text("Quantities:", [ [quant_list[i], np.shape(ave_quant_data[i])] for i in range(len(ave_quant_data))])

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
            "Er":['crimson', r'$P_{rad}$', Er_to_Pr*P_to_cgs, True]
           }

fig, ax1 = plt.subplots(1, 1)

############################# Get log_q for the following lines ########################################
if '_' in quant_list[0]:
    color, label_q, q_fac, log_q = plot_dic[quant_list[0].split('_')[0]]    
else:
    color, label_q, q_fac, log_q = plot_dic[quant_list[0]]

################################# Use log_q to get y_label #############################################
y_label = " "
if log_q:
    y_label += r'$\log_{10}$'+'('
if quant_list[0][0] == 'P' or quant_list[0][0] == 'p' or quant_list[0][0] == 'r' or quant_list[0][0] == 'S':
    y_label += 'Pressure'
if quant_list[0][0] == 'B':
    y_label += 'Magnetic Fields'
if log_q:
    y_label += ')'
if verbose == "False": Print_subtitle("Y-label of the plot:", y_label)
########################################################################################################


if x_axis == "vs_theta":
    t_start = str(int(t[0]*time_to_sec))
    t_end = str(int(t[-1]*time_to_sec))
    ax1.set_title(str(case)+' between t='+t_start+' and t='+t_end+' at r=' + str(radius_select), fontdict = font)
    x_q = th
    x_fac = rad_to_deg
    xlabel = r'$\theta$'
    ylabel = y_label
    plt.xlim(45,135)
    #plt.ylim(-4, 4)
    savename = save_path + 'Time_Averaged_'+str(quant_list[0])+'_'+str(case)+'_r'+radius_select+'_'+t_start+'_'+t_end+'.png' 
elif x_axis == "vs_t":
    ax1.set_title("Shell averaged vs Time" + str(case)+ 'at r=' + str(radius_select), fontdict = font) 
    x_q = t
    x_fac = time_to_sec
    xlabel = r'$t (sec)$'
    ylabel = y_label      
    savename = save_path + 'Shell_Averaged_vs_T'+str(quant_list[0])+'_'+str(case)+'_r'+radius_select+'.png' 


for q_idx, item in enumerate(quant_list):
    if '_' in item:   
        color, label_q, q_fac, log_q = plot_dic[item.split('_')[0]] 
    else:                                                                                                                                                                                                               
        color, label_q, q_fac, log_q = plot_dic[item]
    if verbose == "False": Print_text("Plotting quantity ", item, " using color: ", color, " with label: ", label_q, " unit factor: ", q_fac, " log scale? ", log_q)
    Plotting_1D(ax1, x_q, x_fac, False, ave_quant_data[q_idx], q_fac, log_q, color, xlabel, ylabel, label_q, font)  

if verbose == "False": Print_subtitle("Save the plot as:", savename)    
plt.savefig(savename,format='png',dpi=300)

Print_title("========================== Done plotting_1D_tave_vs_theta.py ==========================")
