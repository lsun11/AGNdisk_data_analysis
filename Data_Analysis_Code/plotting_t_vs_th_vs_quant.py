############################################################################################################################################
#THIS SCRIPT PLOTS TIME(X), THETA(Y), VS QUANTUTY, AT A GIVEN RADIUS
############################################################################################################################################
#Arguments:                                                                                                                                                                                                         
# 1. case (Wedge8, Wedge9B, Wedge10)  2. quantity1  3. quantity2  4. quantity3
# 5. select radius  6. plot photosphere  7. start iteration  8. end iteration           
# 9. photosphere mode (-1: pure Rossland/Planck, 1: using sqrt( (kappa-kappaes)*kappaes), 2: using sqrt(kappa*kappaes))
# 10. read/loading time files? 11. succinct output? 


# Example command1: make plots of density vs tims vs theta at radius= 120 rg with photosphere using sqrt(kappa*kappaes) for photosphere
# python plotting_t_vs_th_vs_quant.py Wedge8_2 rho sigma_p None 120 Y 3000 3500 2 True False

# Example command2: Make a series of plots in command1 at r = 120, 140, 160, 180 and 200 rg:
# for i in `seq 120 20 200`; do python plotting_t_vs_th_vs_quant.py Wedge8_2 rho sigma_p None $i Y 2000 3500 -1 True False; done 

# Example command3: Output ONLY quantity(e.g. rho) checkpoint files: 
# python plotting_t_vs_th_vs_quant.py Wedge8_2 rho None None 140 N 2000 3500 -1 False False

# Example command4: Output ONLY Time checkpoint files: 
# python plotting_t_vs_th_vs_quant.py Wedge8_2 None None None 180 N 2000 3500 -1 True False

# Example command5: make plots of radiation temperature vs time vs theta at radius= 120 rg with photosphere using sigma (Rossland) / sigma_p (Planck)
# python plotting_t_vs_th_vs_quant.py Wedge8 Temp_r sigma None 120 Y 2000 3658 -1 True False 

# Example command6: make plots of gas temperature vs time vs theta at radius= 120 rg with photosphere using sigma (Rossland) / sigma_p (Planck)
# python plotting_t_vs_th_vs_quant.py Wedge8 Temp_g sigma None 120 Y 2000 3658 -1 True False 

# Example command7: make plots of BrBphi - rho_rrho_phi vs time vs theta at radius= 120 rg with photosphere using sigma (Rossland)                                            
# python plotting_t_vs_th_vs_quant.py Wedge8 Ang_rp sigma None 120 Y 2000 3658 -1 True False

# Example command8: make plots of opacity kappa/kappa_p vs times vs theta at radius= 120 rg with photosphere using sigma (Rossland)                         
# python plotting_t_vs_th_vs_quant.py Wedge8 kappa None None 120 Y 2000 3658 -1 True False
# python plotting_t_vs_th_vs_quant.py Wedge8 kappa_p None None 120 Y 2000 3658 -1 True False 


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
                                                                                                                                                                                                                                       
font = {'family': 'sans-serif',                                                                                                                                                                                                        
        'color':  'black',                                                                                                                                                                                                             
        'weight': 'normal',                                                                                                                                                                                                            
        'size': 16,                                                                                                                                                                                                                    
        }                                                                                                                                                                                                                              

Print_title("========================== Starting plotting_t_vs_th_vs_quant.py =============================")    


########################################################################################################            
# input arguments                                                                                                   
########################################################################################################            
Input_Arguments(argv)

case = argv[1]
quant1 = argv[2]
quant2 = argv[3]
quant3 = argv[4]
radius_select = argv[5]
plot_phot = argv[6]     # "Y" or "N"
start = int(argv[7])
end   = int(argv[8])
ph_mode = int(argv[9])
read_time = argv[10]
succinct   = argv[11]

########################################################################################################            
#  set paths based on cases                                                                                         
########################################################################################################
if str(case).find('Wedge8') == 0:                                                                                     
    Print_subtitle('We are analyzing the case: Wedge8') 
    checkpoint_path = Checkpoint_dir  + 'Wedge8/'                                                                     
    save_path = Plot_dir + 'Wedge8/'                                                                              
elif str(case).find('Wedge9B') == 0:                                                                                  
    Print_subtitle('We are analyzing the case: Wedge9B') 
    checkpoint_path = Checkpoint_dir  + 'Wedge9B/'                                                                    
    save_path = Plot_dir + 'Wedge9B/'                                                                             
elif str(case).find('Wedge9Res') == 0:                                                                                
    Print_subtitle('We are analyzing the case: Wedge9Res')  
    checkpoint_path = Checkpoint_dir  + 'Wedge9Res/'                                                                  
    save_path = Plot_dir + 'Wedge9Res/'                                                                           
elif str(case).find('Wedge10') == 0:                                                                                  
    Print_subtitle('We are analyzing the case: Wedge10') 
    checkpoint_path = Checkpoint_dir  + 'Wedge10/'                                                                    
    save_path = Plot_dir + 'Wedge10/' 
dir = Data_dir + 'DATA'+str(case)


########################################################################################################            
# Getting the radius and theta arrays:   
# Any hist file can be used to generate these:                                                                      
########################################################################################################            
filename = "hist_"+str(start).zfill(5)+".npz"                                                                       
if succinct == "False": Print_subtitle("Getting Radius and Theta arrays, using file:", filename) 
data = np.load(dir+'/'+filename)                                                                                    
                                                                                                                    
r = Get_All_1D('radius', data, -1, dir, succinct)   
th = Get_All_1D('theta', data, -1, dir, succinct)                                                                      


########################################################################################################                        
# Check and Load files                                                                                                          
########################################################################################################    
quant_list = []
comp_Tr = 0
comp_Tg = 0  
comp_angmom_rp = 0
comp_angmom_tp = 0
comp_kappa = 0

if quant1 == "Temp_r":
    comp_Tr = 1
    quant1 = "Er"
elif quant2 == "Temp_r":    
    comp_Tr = 1
    quant2 = "Er"
if quant1 == "Temp_g":                                                                                                                                  
    if str(radius_select) not in pg_idx:
        Print_subtitle("Wrong radius for pgas, exit!!!")
        os.sys.exit(0)
    else:
        comp_Tg = 1                                                                                                                                         
        quant1 = "rho"                                                                                                                                       
        quant3 = "pg"
if quant1 == "pg" and str(radius_select) not in pg_idx:
    Print_subtitle("Wrong radius for pgas, exit!!!")                                                                                                
    os.sys.exit(0) 
if quant2 == "Temp_g" or quant3 == "Temp_g":
    Print_subtitle("Please set Temp_g to quant1!!!!")                                                                                         
    os.sys.exit(0)
if quant1 == "Ang_rp":                                                                                            
    comp_angmom_rp = 1                                                                                                  
    quant1 = "rhovrvphi"                                                                                               
    quant3 = "BrBp"
if quant1 == "Ang_tp":                                                                                          
    comp_angmom_tp = 1                                                                                          
    quant1 = "rhovtvphi"                                                                                        
    quant3 = "BtBp"
if quant2 == "Ang_rp" or quant2 == "Ang_tp":
    Print_subtitle("Please set Ang_rp or Ang_tp to quant1!!!!")
    os.sys.exit(0)
if quant1 == "kappa":
    comp_kappa = 1
    quant1 = "rho"
    quant2 = "sigma"
if quant1 == "kappa_p":                                                                                                                                       
    comp_kappa = 1                                                                                                                                          
    quant1 = "rho"                                                                                                                                          
    quant2 = "sigma_p"


if quant1 != "None":
   quant_list.append(quant1)
if quant2 != "None":                                                                                                
    quant_list.append(quant2)
if quant3 != "None":                                                                                                                                                 
    quant_list.append(quant3)
if plot_phot == "Y":
    quant_list.extend(["iphot_upper", "iphot_lower"])

if succinct == "False": Print_subsubtitle("The list of quantities we want to work on is:", quant_list)
########################################################################################################          
# Getting the time and quantities arrays:                                                                             
########################################################################################################
t 	    = []
quant1_data = []
quant2_data = []
quant3_data = []
iphot_upper_total = []
iphot_lower_total = []
file_exist = 0


t_filenames_pre = checkpoint_path + "t_vs_th_vs_quant_time_"+str(radius_select)+"__" 
quant_filenames_pre = checkpoint_path + "t_vs_th_vs_quant_"
idx = 0
if plot_phot == "Y":
    if quant1 != "None":
        if succinct == "False": Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files") 
        filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"    
        quant1_data = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))
        idx += 1
    if quant2 != "None":
        if succinct == "False": Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                             
        filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"             
        quant2_data = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct)) 
        idx += 1
    if quant3 != "None":                                                                                                                                             
        if succinct == "False": Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                                                    
        filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"                                                                           
        quant3_data = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))                                         
        idx += 1
    if succinct == "False": Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                             
    #print(idx, str(quant_list[idx]))
    filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"             
    iphot_upper_total = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))
    idx += 1 
    print(idx, str(quant_list[idx]))
    if succinct == "False": Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                             
    filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"             
    t, iphot_lower_total, file_exist, save_start, start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, True, read_time, succinct)) 

    if succinct == "False": 
        Print_subtitle("Data structure after loading saved files:")                                               
        Print_text("Time array:", np.shape(t)[0])                                                                 
        Print_text("Quant_list:", quant_list)
        Print_text("Quantity and its data set:")
        Print_text(quant_list[0]+":", np.shape(quant1_data), quant_list[1]+":", np.shape(quant2_data))
        #Print_text(quant_list[2]+":", np.shape(iphot_upper_total), quant_list[3]+":", np.shape(iphot_lower_total))
        Print_text("Starting iteration existed in the saved data:", save_start)                                   
        Print_text("(!) Starting iteration for new computation:", start)                                          
        Print_text("(!) Expected to read and plot data until iteration:", end) 
else:    
    if len(quant_list) == 0:  # Probably we just want the time data here
        filenames = "It_doenst_matter"
    else:
        filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"             

    if quant2 != "None":
        Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")
        quant1_data = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))    
        idx += 1 
        Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                             
        filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*" 
        if quant3 != "None":
            quant2_data = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))  
            idx += 1                                                                                                                                                                  
            Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")    
            filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"    
            t, quant3_data, file_exist, save_start, start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, True, read_time, succinct))  
        else:
            t, quant2_data, file_exist, save_start, start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, True, read_time, succinct))
    elif quant3 != "None":
        Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")
        quant1_data = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))                                                      
        idx += 1       
        Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                     
        filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*" 
        t, quant3_data, file_exist, save_start, start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, True, read_time, succinct))
    else:
        t, quant1_data, file_exist, save_start, start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, True, read_time, succinct)) 

    if succinct == "False": 
        Print_subtitle("Data structure after loading saved files:")                                                                                               
        Print_text("Time array:", np.shape(t)[0])                                                                                                                 
        Print_text("Quantity and its data set:")                                                                                                                  
       #Print_text(quant_list[0]+":", np.shape(quant1_data), quant_list[1]+":", np.shape(quant2_data))                                                            
       #Print_text(quant_list[2]+":", np.shape(iphot_upper_total), quant_list[3]+":", np.shape(iphot_lower_total))                                                
        Print_text("Starting iteration existed in the saved data:", save_start)                                                                                   
        Print_text("(!) Starting iteration for new computation:", start)                                                                                          
        Print_text("(!) Expected to read and plot data until iteration:", end) 
                                                                                                                                
save_step = 20                                                                                                                  
save_time_file_pre = checkpoint_path + "t_vs_th_vs_quant_time_" + str(radius_select)+"_" 

Print_title(sp20*5) 
########################################################################################################          
# Computing the new time and quantities arrays:                                                                         
########################################################################################################
Print_subtitle("Start Reading New Data to Arrays!!! Total files to read: ", end-start+1) 
pbar = tqdm(total=end-start+1)
for iter in range(start, end):
    filename = "hist_"+str(iter).zfill(5)+".npz"
    if os.path.exists(dir+'/'+filename):
        data = np.load(dir+'/'+filename)
        if succinct == "False":  
            Print_subtitle("Reading file:", filename) 
        else:
            pbar.update(n=1)
        if file_exist == 1:                                                                                                         
            t = np.append(t, Get_Time(data, dir))                                                                                   
        else:                                                                                                                       
            t.append(Get_Time(data, dir))       

        if radius_select != "None":
            if succinct == "False": Print_subsubtitle("Reading 2D data, need to load values at a specific radius:", radius_select)
            shift_rad = [np.abs(rad - int(radius_select)) for rad in r]        
            radius_idx = shift_rad.index(min(shift_rad))        
            if file_exist == 0:
                if quant1 != "None":
                    if quant1 == "pg":
                        quant1_data.append( data[quant1][pg_idx[str(radius_select)]])
                    else:                                                                                            
                        quant1_data.append( [q[radius_idx] for q in data[quant1]])
                if quant2 != "None":
                    quant2_data.append( [q[radius_idx] for q in data[quant2]])
                if quant3 != "None":                                                                                                                                 
                    if quant3 == "pg":
                        quant3_data.append( data[quant3][pg_idx[str(radius_select)]]) 
                    else:
                        quant3_data.append( [q[radius_idx] for q in data[quant3]])
            else:
                if quant1 != "None":
                    if quant1 == "pg":                                         
                        quant1_data = np.vstack( [quant1_data, [data[quant1][pg_idx[str(radius_select)]]] ] )  
                    else:                                                                                                                           
                        quant1_data = np.vstack( [quant1_data, [q[radius_idx] for q in data[quant1]] ] ) 
                if quant2 != "None":  
                    quant2_data = np.vstack( [quant2_data, [q[radius_idx] for q in data[quant2]] ] )
                if quant3 != "None":                                                                                                                                 
                    if quant3 == "pg":
                        quant3_data = np.vstack( [quant3_data, [data[quant3][pg_idx[str(radius_select)]]] ] ) 
                    else:
                        quant3_data = np.vstack( [quant3_data, [q[radius_idx] for q in data[quant3]] ] )
            if succinct == "False": Print_text("Data structure of quantities:", quant1+":", np.shape(quant1_data), quant2+":", np.shape(quant2_data), quant3+":", np.shape(quant3_data))

            if plot_phot == "Y":     # make sure your data[quant2] is kappa (or kappa*rho in reality) here!!!!
                if succinct == "False": Print_text("Finding the photosphere!")
                if ph_mode < 0:
                    iphot_upper = Get_Photosphere(th, [q[radius_idx] for q in data[quant2]], float(radius_select), 'upper')
                    iphot_lower = Get_Photosphere(th, [q[radius_idx] for q in data[quant2]], float(radius_select), 'lower') 
                else:
                    iphot_upper = Get_Photosphere2(th, [q[radius_idx] for q in data[quant2]], [q[radius_idx] for q in data[quant1]], float(radius_select), ph_mode, 'upper')
                    iphot_lower = Get_Photosphere2(th, [q[radius_idx] for q in data[quant2]], [q[radius_idx] for q in data[quant1]], float(radius_select), ph_mode, 'lower')
                if file_exist == 0:
                    iphot_upper_total.append(iphot_upper)
                    iphot_lower_total.append(iphot_lower)
                else:
                    iphot_upper_total = np.append(iphot_upper_total, iphot_upper)
                    iphot_lower_total = np.append(iphot_lower_total, iphot_lower)
            if succinct == "False": Print_text("Data structure of quantities:","iphot_upper_total:",np.shape(iphot_upper_total), "iphot_lower_total:",np.shape(iphot_lower_total))
        else:
            if succinct == "False": Print_subsubtitle("Reading 1D data, radius already specified")
            if file_exist == 0:
                if quant1 != "None":  
                    quant1_data.append(data[quant1])
                if quant2 != "None":  
                    quant2_data.append(data[quant2])
                if quant3 != "None":                                                                                                                                 
                    quant3_data.append(data[quant3])
            else:
                if quant1 != "None":  
                    quant1_data = np.vstack( [quant1_data, data[quant1]])
                if quant2 != "None":   
                    quant2_data = np.vstack( [quant2_data, data[quant2]])
                if quant3 != "None":                                                                                                                                 
                    quant3_data = np.vstack( [quant2_data, data[quant3]])
            if succinct == "False": Print_text("Data structure of quantities:", quant1+":", np.shape(quant1_data), quant2+":",np.shape(quant2_data), quant3+":", np.shape(quant3_data)) 
            if plot_phot == "Y":     # make sure your data[quant2] is kappa (or kappa*rho in reality) here!!!!                                                                                            
                Print_text("Finding the photosphere!")  
                if ph_mode < 0: 
                    iphot_upper = Get_Photosphere(th, data[quant2], float(radius_select), 'upper')
                    iphot_lower = Get_Photosphere(th, data[quant2], float(radius_select), 'lower')
                else:
                    iphot_upper = Get_Photosphere(th, data[quant2], float(radius_select), ph_mode, 'upper')                                                                                          
                    iphot_lower = Get_Photosphere(th, data[quant2], float(radius_select), ph_mode, 'lower')
                iphot_upper_total.append(iphot_upper)                                                                                                                     
                iphot_lower_total.append(iphot_lower) 
            if succinct == "False": Print_text("Data structure of quantities:","iphot_upper_total:",np.shape(iphot_upper_total), "iphot_lower_total:",np.shape(iphot_lower_total))

########################################################################################################                      
# Save files since appending above is too slow!!!                                                                             
# TODO: move save_step out of Save_Files
########################################################################################################
        if iter > start and (iter%save_step == 0 or iter == end-1):  
            idx = 0
            if quant1 != "None":  
                if succinct == "False": Print_subtitle("Saving files At Iteration", iter, "quantity: ", quant_list[idx]) 
                Save_Files(save_step, iter, save_start, start, end, quant1_data, quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_", succinct)
                idx += 1
            if quant2 != "None":
                if succinct == "False": Print_subtitle("Saving files At Iteration", iter, "quantity: ", quant_list[idx])                                                                                          
                Save_Files(save_step, iter, save_start, start, end, quant2_data, quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_", succinct) 
                idx += 1
            if quant3 != "None":                                                                                                                                     
                if succinct == "False": Print_subtitle("Saving files At Iteration", iter, "quantity: ", quant_list[idx])                                             
                Save_Files(save_step, iter, save_start, start, end, quant3_data, quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_", succinct)    
                idx += 1
            if plot_phot == "Y":
                if succinct == "False": Print_subtitle("Saving files At Iteration", iter, "quantity: ", quant_list[idx])                                                                                          
                Save_Files(save_step, iter, save_start, start, end, iphot_upper_total, quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_", succinct)
                idx += 1 
                if succinct == "False": Print_subtitle("Saving files At Iteration", iter, "quantity: ", quant_list[idx])                                                                                  
                Save_Files(save_step, iter, save_start, start, end, iphot_lower_total, quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_", succinct)
            if read_time == "True":
                if succinct == "False": Print_subtitle("Saving time At Iteration", iter)
                Save_Files(save_step, iter, save_start, start, end, t, save_time_file_pre, succinct)

if succinct == "False": Print_subtitle("Done reading data!")                                                                                
if succinct == "False": Print_text("The shape of the datasets are:", np.shape(quant1_data), np.shape(quant2_data), np.shape(quant3_data))



########################################################################################################
# Convert Er to Temperature using Er = aT^4  
########################################################################################################
if comp_Tr == 1:
    if succinct == "False": Print_subtitle("Computing Radiation Temperature!!!")
    if quant1 == 'Er':
        quant1_data = [[qq**(0.25) for qq in q] for q in quant1_data]
    elif quant2 == 'Er':                                                                                                   
        quant2_data = [[qq**(0.25) for qq in q] for q in quant2_data] 
    quant1 = "Temp_r"

if comp_Tg == 1:                                                                                                                                        
    if succinct == "False": Print_subtitle("Computing Gas Temperature!!!") 
    quant1_data = [[qq3/qq1 for qq1,qq3 in zip(q1, q3)] for q1, q3 in zip(quant1_data, quant3_data)]    
    quant1 = "Temp_g" 
    
if comp_angmom_rp == 1:
   if succinct == "False": Print_subtitle("Computing rhovrvphi - BrBp!!!") 
   print(np.shape(quant1_data), np.shape(quant3_data))
   quant1_data = [[qq1 - qq3 for qq1,qq3 in zip(q1, q3)] for q1, q3 in zip(quant1_data, quant3_data)]
   quant1 = "Ang_rp"
if comp_angmom_tp == 1:                                                                                         
   if succinct == "False": Print_subtitle("Computing rhovtvphi - BtBp!!!")                                      
   quant1_data = [[qq1 - qq3 for qq1,qq3 in zip(q1, q3)] for q1, q3 in zip(quant1_data, quant3_data)] 
   quant1 = "Ang_tp" 

if comp_kappa == 1:
   if succinct == "False": Print_subtitle("Computing Opacity")
   # Note: the quantity in the line below: sigma(code)/rho(code) is just kappa/kappa_es! No additional factor needed!
   quant1_data = [[qq2/qq1 for qq1,qq2 in zip(q1, q2)] for q1, q2 in zip(quant1_data, quant2_data)] 
   #quant1_data = [[qq2/kappaes_code for qq2 in q2] for q2 in quant2_data] 

Print_title(sp30, "Start Plotting!!!!", sp30)                                                                     
########################################################################################################          
# Plotting and save              
########################################################################################################  
v_max = None
v_min = None


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
    if (str(quant2)[0:5]) == "sigma":             
        color = 'bone' 
        if (str(quant2)[-1]) == "p":
            #pre_str = '\\kappa_{P}/\\kappa_{es}'
            #cbar_str = r'$\kappa_{Planck}/\kappa_{es}$'
            pre_str = '\\kappa_{Plank}'
            cbar_str = r'$\kappa_{Planck}$'
            quant1 = "kappa_p"
        else:
            #pre_str = '\\kappa_{R}/\\kappa_{es}' 
            #cbar_str = r'$\kappa_{Rossland}/\kappa_{es}$' 
            pre_str = '\\kappa_{Rossland}'                                                                                                                                                                             
            cbar_str = r'$\kappa_{Rossland}$'
            quant1 = "kappa"
        fac  = 1.0/(dens_to_cgs*r_to_cgs)
        #v_max = 3.0
        #v_min = 0.0
        #fac = 1.0
        logscale = False
    else:
        color = 'BrBG_r'                                                                                             
        pre_str = 'log_{10}(\\rho/\\rho_0)'
        cbar_str = r'$log_{10}(\rho/\rho_0)$'                                                                                   
        #fac  = 1.0/dens_to_cgs                                                                                               
        fac  = 1.0          # code unit rho is just rho_cgs/rho_0
        logscale = True                                                                                              
    phot_color = 'yellow'
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
elif str(str(quant1)[0:6]) == "Temp_r": # no _2          
    color = 'Spectral'                                    
    pre_str = 'T_{rad} (K)'                                                                                                  
    cbar_str = r'$\log_{10}(T_{rad}) (K)$'                                                                                               
    fac = (P_to_cgs/rad_const)**0.25
    logscale = True
    phot_color = 'black'
elif str(str(quant1)[0:6]) == "Temp_g": # no _2                                                                                                         
    color = 'pink'                                                                                                                                  
    pre_str = 'T_{gas} (K)'                                                                                                                                   
    cbar_str = r'$\log_{10}(T_{gas}) (K)$'                                                                                                                    
    fac = 1.0e5
    logscale = True                                                                                                                                     
    phot_color = 'black' 
    v_max = 6.0
    v_min = 2.0
elif str(str(quant1)[0:6]) == "Ang_rp":
    color = 'winter'
    pre_str = 'log_{10}(\\rho_r\\rho_{\\phi} - B_rB_{\\phi})'
    cbar_str = r'$\log_{10}(\rho_r\rho_{\phi} - B_rB_{\phi})$'
    fac = P_to_cgs
    logscale = True 
    phot_color = 'white' 
elif str(str(quant1)[0:6]) == "Ang_tp":
    color = 'cool'                                                                                            
    pre_str = 'log_{10}(\\rho_{\\theta}\\rho_{\\phi} - B_{\\theta}B_{\\phi})'   
    cbar_str = r'$\log_{10}(\rho_{\theta}\rho_{\phi} - B_{\theta}B_{\phi})$'                                                             
    fac = P_to_cgs                                                                                       
    logscale = True                                                                                            
    phot_color = 'white'
print (str(quant1))


Theta, Time = np.meshgrid(th,t)

Theta *= 57.2958 
Time *= 495.0766                   #convert to second                                                                                                                                                                                                                                                                                                                                                                                                             
t_phys =[time * 495.0766 for time in t]    
start_time = int(t_phys[0])    
end_time = int(t_phys[-1]) 


radius = {'1': 120, '2':140, '3':160, '4':180}                                                                                                                                                                                   
if str(quant1)[-1]  in radius.keys():                                                                                                                                                                                                        
    string_radius = str(radius[str(quant1)[-1]])                                                                                                                                                                                            
else:                                                                                                                                                                                                                            
    string_radius = str(int(radius_select))
title = r'${}$ at $r =$'.format(pre_str) + string_radius +'$ r_g$' + ' from ' + str(start_time) + ' to ' + str(end_time)
xlabel = r'$t (s)$'
ylabel = r'$\theta(^o)$'



################## Logscale of q has been taken care of, set log_q to False here ###################################################################
ax1 = Plotting_Mesh_YX(t, th, np.abs(quant1_data), time_to_sec, rad_to_deg, fac, False, False, logscale, color, None, None, v_max, v_min, title, xlabel, ylabel, cbar_str, font)

savename = save_path + str(quant1) + '_vs_t_and_theta_' + str(case) + '_' + string_radius + '_' + str(save_start) + '_' +str(end) + '_' + 'plot.png' 
################## Plot the photosphere curve ######################################################################################################
if plot_phot == "Y":
    phot_th_upper = [th[it]*rad_to_deg for it in iphot_upper_total]                                                                                                                                                                  
    phot_th_lower = [th[it]*rad_to_deg for it in iphot_lower_total]
    quant_list = [phot_th_upper, phot_th_lower]
    Plotting_1D(ax1, t, time_to_sec, False, quant_list, 1.0, False, phot_color, None, None, None, font)
    ### Add an extra part of the filename
    split_name = savename.split('_')
    split_name.insert(-1, 'photosphere')
    savename = ('_').join(split_name)
################## Save the plot ###################################################################################################################
plt.savefig(savename,format='png',dpi=300)
#####################################################################################################################################################


Print_title("========================== Done plotting_t_vs_th_vs_quant.py =============================")
