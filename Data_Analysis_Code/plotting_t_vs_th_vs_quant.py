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
quant_list_data  = []
save_start_list  = []
save_end_list    = []
read_start_list  = []
read_end_list    = []
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
        quant1_data, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))
        quant_list_data.append(quant1_data)
        save_start_list.append(save_start)
        read_start_list.append(read_start)
        idx += 1
    if quant2 != "None":
        if succinct == "False": Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                             
        filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"             
        quant2_data, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct)) 
        quant_list_data.append(quant2_data)
        save_start_list.append(save_start)                                                                                                                    
        read_start_list.append(read_start)
        idx += 1
    if quant3 != "None":                                                                                                                                             
        if succinct == "False": Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                                                    
        filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"                                                                           
        quant3_data, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))                                         
        quant_list_data.append(quant3_data)
        save_start_list.append(save_start)                                                                                                                    
        read_start_list.append(read_start)
        idx += 1
    if succinct == "False": Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                             
    filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"             
    iphot_upper_total, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))
    quant_list_data.append(iphot_upper_total)
    save_start_list.append(save_start)                                                                                                                    
    read_start_list.append(read_start)
    idx += 1 

    if succinct == "False": Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                             
    filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"             
    #t, iphot_lower_total, file_exist, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, True, read_time, succinct)) 
    iphot_lower_total, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))
    quant_list_data.append(iphot_lower_total) 
    save_start_list.append(save_start)                                                                                                                        
    read_start_list.append(read_start)

    if read_time == "True":
        if succinct == "False": Print_subtitle("Compute time!!! First check saved files")
        t_filenames = t_filenames_pre + "*"
        t, time_save_start, time_read_start = list(Check_Load_Files(t_filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))
        quant_list_data.append(t)



    start = np.min(save_start_list)
    if succinct == "False": 
        Print_subtitle("Data structure after loading saved files:")                                               
        Print_text("Time array length:", np.shape(t)[0]) 
        Print_text("The list of quantities and the shape of their data set:")
        Print_text( [ (q, np.shape(quant_list_data[quant_list.index(q)])) for q in quant_list] )
        Print_text("Starting iteration existed in the saved data:", save_start)                                   
        Print_text("Save start and read start lists are:", save_start_list, read_start_list)
        Print_text("(!) Starting iteration for new computation:", np.min(read_start_list))                                          
        Print_text("(!) Expected to read and plot data until iteration:", end) 

else:    
    if len(quant_list) == 0:  # Probably we just want the time data here
        filenames = "It_doenst_matter"
    else:
        filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"             

    if quant2 != "None":
        Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")
        quant1_data, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))    
        quant_list_data.append(quant1_data)
        save_start_list.append(save_start)                                                                                                                        
        read_start_list.append(read_start) 
        idx += 1 
        Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                             
        filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*" 
        if quant3 != "None":
            quant2_data, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))  
            quant_list_data.append(quant2_data)
            save_start_list.append(save_start)
            read_start_list.append(read_start)                                                                                                                         
            idx += 1                                                                                                                                                                  
            Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")    
            filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*"    
            quant3_data, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))  
            quant_list_data.append(quant3_data)
            save_start_list.append(save_start)                                                                                                                        
            read_start_list.append(read_start) 
        else:
            quant2_data, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))
            quant_list_data.append(quant2_data)
            save_start_list.append(save_start)                                                                                                                        
            read_start_list.append(read_start) 
    elif quant3 != "None":
        Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")
        quant1_data, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))                                                      
        quant_list_data.append(quant1_data)
        save_start_list.append(save_start)
        read_start_list.append(read_start)                                                                                                                         
        idx += 1       
        Print_subtitle("Compute " + str(quant_list[idx]) + "!!! First check saved files")                     
        filenames = quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_*" 
        quant3_data, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))
        quant_list_data.append(quant3_data)
        save_start_list.append(save_start)                                                                                                                        
        read_start_list.append(read_start) 
    else:
        quant_list_data.append(quant1_data)
        quant1_data, save_start, read_start = list(Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct)) 
        save_start_list.append(save_start)
        read_start_list.append(read_start)                                                                                                                         
   
    if read_time == "True":                                                                                                                                                                                                                             
        if succinct == "False": Print_subtitle("Compute time!!! First check saved files") 
        t_filenames = t_filenames_pre + "*"
        t, time_save_start, time_read_start = list(Check_Load_Files(t_filenames, t_filenames_pre, file_exist, start, end, False, read_time, succinct))   
        quant_list_data.append(t)

    start = np.min(save_start_list)
    if succinct == "False": 
        Print_subtitle("Data structure after loading saved files:")                                                                                               
        Print_text("Time array length:", np.shape(t)[0])                                                              
        Print_text("The list of quantities and the shape of their data set:")
        Print_text( [ (q, np.shape(quant_list_data[quant_list.index(q)])) for q in quant_list] ) 
        Print_text("Starting iteration existed in the saved data:", save_start)                                          
        Print_text("Save start and read start lists are:", save_start_list, read_start_list)                             
        Print_text("(!) Starting iteration for new computation:", start)
        Print_text("(!) Expected to read and plot data until iteration:", end)

                    

save_step = 20                                                                                                                  
save_time_file_pre = checkpoint_path + "t_vs_th_vs_quant_time_" + str(radius_select)+"_" 

Print_title(sp20*5) 
########################################################################################################          
# Computing the new time and quantities arrays:                                                                         
########################################################################################################
Print_subtitle("Start Reading New Data to Arrays!!! Total files to read: ", end-start+1) 
pbar = tqdm(total=end-start)
for iter in range(start, end):
    filename = "hist_"+str(iter).zfill(5)+".npz"
    if os.path.exists(dir+'/'+filename):
        data = np.load(dir+'/'+filename)
        if succinct == "False":  
            Print_subtitle("Reading file:", filename) 
        else:
            pbar.update(n=1)


        if iter >= time_read_start:
            if len(t) > 0:   
                t = np.append(t, Get_Time(data, dir))                                                                                   
            else:                                                                                                                       
                t.append(Get_Time(data, dir))       

        
        if radius_select != "None":
            if succinct == "False": Print_subsubtitle("Reading 2D data, need to load values at a specific radius:", radius_select)
            shift_rad = [np.abs(rad - int(radius_select)) for rad in r]        
            radius_idx = shift_rad.index(min(shift_rad))        

            if quant1 != "None" and iter >= read_start_list[quant_list.index(quant1)]:
                Print_text("Loading quantity: ", str(quant1))
                if np.shape( quant_list_data[quant_list.index(quant1)] )[0] == 0:
                    if quant1 == "pg":
                        quant_list_data[quant_list.index(quant1)].append( data[quant1][pg_idx[str(radius_select)]] )
                    else:                                                                                            
                        quant_list_data[quant_list.index(quant1)].append( [q[radius_idx] for q in data[quant1]] )
                else:
                    if quant1 == "pg":                                                                                                                                                               
                        quant_list_data[quant_list.index(quant1)] =  np.vstack( [quant_list_data[quant_list.index(quant1)], [data[quant1][pg_idx[str(radius_select)]]] ] )        
                    else:                                                                                                                                                                                 
                        quant_list_data[quant_list.index(quant1)] =  np.vstack( [quant_list_data[quant_list.index(quant1)], [q[radius_idx] for q in data[quant1]] ] )  
 
            if quant2 != "None" and iter >= read_start_list[quant_list.index(quant2)]:
                Print_text("Loading quantity: ", str(quant2))  
                if np.shape( quant_list_data[quant_list.index(quant2)] )[0] == 0: 
                    quant_list_data[quant_list.index(quant2)].append( [q[radius_idx] for q in data[quant1]] )
                else:
                    quant_list_data[quant_list.index(quant2)] =  np.vstack( [quant_list_data[quant_list.index(quant2)], [q[radius_idx] for q in data[quant2]] ] )

            if quant3 != "None" and iter >= read_start_list[quant_list.index(quant3)]:
                Print_text("Loading quantity: ", str(quant3))
                if np.shape( quant_list_data[quant_list.index(quant3)] )[0] == 0: 
                    if quant3 == "pg":
                        quant_list_data[quant_list.index(quant3)].append( data[quant3][pg_idx[str(radius_select)]] ) 
                    else:
                        quant_list_data[quant_list.index(quant3)].append( [q[radius_idx] for q in data[quant3]] )
                else:
                    if quant3 == "pg":
                        quant_list_data[quant_list.index(quant3)] =  np.vstack( [quant_list_data[quant_list.index(quant3)], [data[quant3][pg_idx[str(radius_select)]]] ] ) 
                    else:
                        quant_list_data[quant_list.index(quant3)] =  np.vstack( [quant_list_data[quant_list.index(quant3)], [q[radius_idx] for q in data[quant3]] ] )
            if plot_phot == "Y" and iter >= read_start_list[quant_list.index("iphot_upper")]:
                if succinct == "False": Print_text("Finding the photosphere!")
                if ph_mode < 0:
                    iphot_upper = Get_Photosphere(th, [q[radius_idx] for q in data[quant2]], float(radius_select), 'upper')
                    iphot_lower = Get_Photosphere(th, [q[radius_idx] for q in data[quant2]], float(radius_select), 'lower') 
                else:
                    iphot_upper = Get_Photosphere2(th, [q[radius_idx] for q in data[quant2]], [q[radius_idx] for q in data[quant1]], float(radius_select), ph_mode, 'upper')
                    iphot_lower = Get_Photosphere2(th, [q[radius_idx] for q in data[quant2]], [q[radius_idx] for q in data[quant1]], float(radius_select), ph_mode, 'lower')
                if  np.shape( quant_list_data[quant_list.index("iphot_upper")] )[0] == 0:
                    quant_list_data[quant_list.index("iphot_upper")].append(iphot_upper)
                else:
                    quant_list_data[quant_list.index("iphot_upper")] = np.append( quant_list_data[quant_list.index("iphot_upper")], iphot_upper)  

                if  np.shape( quant_list_data[quant_list.index("iphot_lower")] )[0] == 0: 
                    quant_list_data[quant_list.index("iphot_lower")].append(iphot_lower)
                else:
                    quant_list_data[quant_list.index("iphot_lower")] = np.append( quant_list_data[quant_list.index("iphot_lower")], iphot_lower)

            if succinct == "False": 
                Print_text("Data structure of quantities:", "time,", len(t), [ (q, np.shape(quant_list_data[quant_list.index(q)]) ) for q in quant_list ] )

        else:  # radius_select == None (probably will be deprecated soon)
            if succinct == "False": Print_subsubtitle("Reading 1D data, radius already specified")
            if file_exist == 0:
                if quant1 != "None" and iter >= read_start_list[quant_list.index(quant1)]:
                    Print_text("Loading quantity: ", str(quant1))  
                    quant_list_data[quant_list.index(quant1)].append( data[quant1] )
                if quant2 != "None" and iter >= read_start_list[quant_list.index(quant2)]:
                    Print_text("Loading quantity: ", str(quant2))  
                    quant_list_data[quant_list.index(quant2)].append( data[quant2] )
                if quant3 != "None" and iter >= read_start_list[quant_list.index(quant3)]:
                    Print_text("Loading quantity: ", str(quant3))                                                                                                                  
                    quant_list_data[quant_list.index(quant3)].append( data[quant3] ) 
            else:
                if quant1 != "None" and iter >= read_start_list[quant_list.index(quant1)]:
                    Print_text("Loading quantity: ", str(quant1))    
                    quant_list_data[quant_list.index(quant1)] = np.vstack( [quant_list_data[quant_list.index(quant1)], data[quant1]] )
                if quant2 != "None" and iter >= read_start_list[quant_list.index(quant2)]:
                    Print_text("Loading quantity: ", str(quant2))  
                    quant_list_data[quant_list.index(quant2)] = np.vstack( [quant_list_data[quant_list.index(quant2)], data[quant2]] ) 
                if quant3 != "None" and iter >= read_start_list[quant_list.index(quant3)]:
                    Print_text("Loading quantity: ", str(quant3))                                                                                                                                   
                    quant_list_data[quant_list.index(quant3)] = np.vstack( [quant_list_data[quant_list.index(quant3)], data[quant3]] )
            if plot_phot == "Y" and iter > read_start_list[quant_list.index("iphot_upper")]:
                Print_text("Finding the photosphere!")  
                if ph_mode < 0: 
                    iphot_upper = Get_Photosphere(th, data[quant2], float(radius_select), 'upper')
                    iphot_lower = Get_Photosphere(th, data[quant2], float(radius_select), 'lower')
                else:
                    iphot_upper = Get_Photosphere2(th, data[quant2], float(radius_select), ph_mode, 'upper')                                                                                          
                    iphot_lower = Get_Photosphere2(th, data[quant2], float(radius_select), ph_mode, 'lower')
                quant_list_data[quant_list.index("iphot_upper")].append(iphot_upper)                                                                                                                                                             
                quant_list_data[quant_list.index("iphot_lower")].append(iphot_lower)
            if succinct == "False": Print_text("Data structure of quantities:", [ (q, np.shape(quant_list_data[quant_list.index(q)])) for q in quant_list])

########################################################################################################                      
# Save files since appending above is too slow!!!                                                                             
# TODO: move save_step out of Save_Files
########################################################################################################
        if iter > start and (iter%save_step == 0 or iter == end-1):  
            idx = 0
            if quant1 != "None":  
                if succinct == "False": Print_subtitle("Saving files At Iteration", iter, "quantity: ", quant_list[idx]) 
                Save_Files(save_step, iter, save_start, start, end, quant_list_data[quant_list.index(quant1)], quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_", succinct)
                idx += 1
            if quant2 != "None":
                if succinct == "False": Print_subtitle("Saving files At Iteration", iter, "quantity: ", quant_list[idx])                                                                                          
                Save_Files(save_step, iter, save_start, start, end, quant_list_data[quant_list.index(quant2)], quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_", succinct) 
                idx += 1
            if quant3 != "None":                                                                                                                                     
                if succinct == "False": Print_subtitle("Saving files At Iteration", iter, "quantity: ", quant_list[idx])                                             
                Save_Files(save_step, iter, save_start, start, end, quant_list_data[quant_list.index(quant3)], quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_", succinct)    
                idx += 1
            if plot_phot == "Y":
                if succinct == "False": Print_subtitle("Saving files At Iteration", iter, "quantity: ", quant_list[idx])                                                                                          
                Save_Files(save_step, iter, save_start, start, end, quant_list_data[quant_list.index("iphot_upper")], quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_", succinct)
                idx += 1 
                if succinct == "False": Print_subtitle("Saving files At Iteration", iter, "quantity: ", quant_list[idx])                                                                                  
                Save_Files(save_step, iter, save_start, start, end, quant_list_data[quant_list.index("iphot_lower")], quant_filenames_pre + str(quant_list[idx])+"_"+str(radius_select)+"_", succinct)
            if read_time == "True":
                if succinct == "False": Print_subtitle("Saving files At Iteration", iter, "quantity: time")
                Save_Files(save_step, iter, save_start, start, end, t, save_time_file_pre, succinct)

if succinct == "False": Print_subtitle("Done reading data!")                                                                                
if succinct == "False": Print_text("The shape of the datasets are:", [ (q, np.shape(quant_list_data[quant_list.index(q)])) for q in quant_list] ) 



########################################################################################################
# Convert Er to Temperature using Er = aT^4  
########################################################################################################
if comp_Tr == 1:
    if succinct == "False": Print_subtitle("Computing Radiation Temperature!!!")
    if quant1 == 'Er':
        quant1_data = [[qq**(0.25) for qq in q] for q in quant_list_data[quant_list.index(quant1)]]
    elif quant2 == 'Er':                                                                                                   
        quant2_data = [[qq**(0.25) for qq in q] for q in quant_list_data[quant_list.index(quant2)]] 
    quant1 = "Temp_r"

if comp_Tg == 1:                                                                                                                                        
    if succinct == "False": Print_subtitle("Computing Gas Temperature!!!") 
    quant1_data = [[qq3/qq1 for qq1,qq3 in zip(q1, q3)] for q1, q3 in zip(quant_list_data[quant_list.index(quant1)], quant_list_data[quant_list.index(quant3)])]    
    quant1 = "Temp_g" 
    
if comp_angmom_rp == 1:
   if succinct == "False": Print_subtitle("Computing rhovrvphi - BrBp!!!") 
   print(np.shape(quant1_data), np.shape(quant3_data))
   quant1_data = [[qq1 - qq3 for qq1,qq3 in zip(q1, q3)] for q1, q3 in zip(quant_list_data[quant_list.index(quant1)], quant_list_data[quant_list.index(quant3)])]
   quant1 = "Ang_rp"
if comp_angmom_tp == 1:                                                                                         
   if succinct == "False": Print_subtitle("Computing rhovtvphi - BtBp!!!")                                      
   quant1_data = [[qq1 - qq3 for qq1,qq3 in zip(q1, q3)] for q1, q3 in zip(quant_list_data[quant_list.index(quant1)], quant_list_data[quant_list.index(quant3)])] 
   quant1 = "Ang_tp" 

if comp_kappa == 1:
   if succinct == "False": Print_subtitle("Computing Opacity")
   # Note: the quantity in the line below: sigma(code)/rho(code) is just kappa/kappa_es! No additional factor needed!
   quant1_data = [[qq2/qq1 for qq1,qq2 in zip(q1, q2)] for q1, q2 in zip(quant_list_data[quant_list.index(quant1)], quant_list_data[quant_list.index(quant2)])] 

if comp_Tr+comp_Tg+comp_angmom_rp+comp_angmom_tp+comp_kappa == 0:
    quant1_data= quant_list_data[quant_list.index(quant1)] 

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
