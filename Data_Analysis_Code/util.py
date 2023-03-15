# This script has all the basic function such as getting a certain quantity vs r or theta, average, etc.

import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import sys

from parameters import *

sp1 = " "
sp2 = "  "
sp3 = "   "
sp20 = sp2*10
sp30 = sp3*10

def print_e(*arg):
    print(*arg,end="")

########################################################################################################                                                                                                            
# Fucntion: Print_title  -- print noticing texts (title)
########################################################################################################                                                                                                            
# Parameters:                                                                                                                                                                                                        
########################################################################################################                                                                                                            
# quant: *arg (whatever the print function passes)
########################################################################################################   
def Print_title(*arg):
    print("\n")  
    boarder = "#" * (len(str(arg[:]))+60)
    print(boarder)
    print(sp30, *arg, sp30)
    print(boarder)
    #print("\n")  


########################################################################################################                                                                                                            
# Fucntion: Print_subtitle  -- print noticing texts (subtitle)                                                                                                                                                            
########################################################################################################                                                                                                            
# Parameters:                                                                                                                                                                                                       
########################################################################################################                                                                                                            
# quant: *arg (whatever the print function passes)                                                                                                                                                                  
########################################################################################################                                                                                                            
def Print_subtitle(*arg):                                                                                                                                                                                              
    print("\n")
    boarder = sp2 + "=" * len(str(arg[:]))      
    print(boarder)                                                                                                                                                                                                  
    print(sp3 + "o", *arg)                                                                                                                                                                                                
    print(boarder)                                                                                                                                                                                                  
    #print("\n") 


def Print_subtitle_e(*arg):                                                                                                
    #print("\n")                                                                                                          
    #boarder = sp2 + "=" * len(str(arg[:]))                                                                               
    #print(boarder)                                                                                                       
    print_e(sp3 + "o", *arg)                                                                                               
    #print(boarder)                                                                                                       
    #print("\n") 

########################################################################################################                                                                                                            
# Fucntion: Print_subsubtitle  -- print noticing texts (subsubtitle)                                                                                                                                                      
########################################################################################################                                                                                                            
# Parameters:                                                                                                                                                                                                       
########################################################################################################                                                                                                            
# quant: *arg (whatever the print function passes)                                                                                                                                                                  
########################################################################################################                                                                                                            
def Print_subsubtitle(*arg):                                                                                                                                                                                           
    boarder = sp2 + "-" * len(str(arg[:]))                                                                                                                                                                                
    print(boarder)   
    print(sp3 + "-", *arg)
    print(boarder)
    #print("\n")                                                                                                                                                                                                    


########################################################################################################                                                                                                            
# Fucntion: Print_text  -- print noticing texts (textsize)        
########################################################################################################                                                                                                            
# Parameters:                                                                                                                                                                                                       
########################################################################################################                                                                                                            
# quant: *arg (whatever the print function passes)                                                                                                                                                                  
########################################################################################################                                                                                                            
def Print_text(*arg):
    print(sp3 + sp2 + "--", *arg) 



########################################################################################################                                             
# Fucntion: Input_Arguments  -- Reading the input arguments
########################################################################################################                                             
# Parameters:                                                                                                                                        
########################################################################################################                                             
# quant: argv (object in module sys to take all the input arguments) 
######################################################################################################## 

def Input_Arguments(Argv):
    Argv_final = []
    for arg in Argv[:]:                                                                                            
        if arg.endswith('.txt'):                                                                                        
           with open(arg[:]) as f:                                                                                     
               file_args = f.read().splitlines()                                                                       
               Argv_final[:] = [arg for arg in file_args if arg[0] != "#"]                                                   
        else:                                                                                               
           Argv_final.append(arg)                                                                                            
                                                                                                                   
    Print_subtitle("Function Input_Arguments: Your input commands are:")
    Print_text(Argv_final)
    return Argv_final



########################################################################################################
# Fucntion: Get_All_1D  -- Get the entire 1D array of radius or theta
######################################################################################################## 
# Parameters:
########################################################################################################
# quant: quantity to get (radius/theta)
# data: the dataset to load the data, generated outside the function 
# cutoff: required if only need part of the data, set it to -1 if we want the entire array
# dir: directory of the data file
# succinct: If we print succinct ('True') info or detailed ('False') info.
########################################################################################################

def Get_All_1D(quant, data, cutoff, dir, succinct):
    result = []
    if succinct == "False": Print_text("Function Get_All_1D: Getting the array of", quant)
    if cutoff > 0:
        result.append(data[quant][:cutoff])
    else:
        result.append(data[quant])
    return result[0]



########################################################################################################
# Fucntion: Get_Time  -- Get the time array for each file
########################################################################################################                 
# Parameters:                             
########################################################################################################
# data: the data set, this function reads external data set to avoid repeated loading
# dir: directory of the data file 
######################################################################################################## 

def Get_Time(data, dir):
    result = []                                                                                                          
    t_data = data['time']
    if len(np.shape(t_data)) == 1:
        result.append(t_data[0])
    else:
        result.append(t_data)
    return result[0]







########################################################################################################                 
# Fucntion: Get_Photosphere  -- Get the theta index/value of the photosphere, at a certain radius                                                                
########################################################################################################                 
# Parameters:                                                                                                            
########################################################################################################
# theta_data: The thetea array, should be from 0 to pi
# kappa_data: The density*opacity array (same size as theta)
# radius:  The radius
# hemi:    Which hemisphere we are integrating (upper or lower)
######################################################################################################## 

def Get_Photosphere(theta_data, kappa_data, radius, hemi):
    #print("Finding theta of photosphere on the " + str(hemi) + " hemisphere!")
    data_len = len(theta_data)
    
    if hemi == 'lower':
        kappa_data = kappa_data[::-1]
        theta_data = theta_data[::-1]

    iphot = 0
    tau = np.zeros(data_len//2)
    for it in range(1, data_len//2, 1):
        #print(it, it-step, tau[it], tau[it-step], radius, kappa_data[it])
        #tau[it] = tau[it-1] + 0.5*radius*kappaes_code*(kappa_data[it-1] + kappa_data[it])*(theta_data[it] - theta_data[it-1]) 
        tau[it] = tau[it-1] + 0.5*radius*1.0*(kappa_data[it-1] + kappa_data[it])*(theta_data[it] - theta_data[it-1])
        if np.abs(tau[it]) >= 1:    # found photosphere, no need to go further                            
           #print("Hit!!!", it, tau[it], tau[it-1], kappa_data[it-1], kappa_data[it])
           iphot = it if (tau[it]+tau[it-1]) < 2 else it-1    
           break                                                                                                         
        iphot = it 
 
    if hemi == 'lower':
        iphot = data_len - iphot -1
    return iphot



######################################################################################################## 
# Fucntion: Get_Photosphere2  -- Get the theta index of the photosphere, at a certain radius (diff kappa) 
######################################################################################################## 
# Parameters: 
######################################################################################################## 
# theta_data: The thetea array, should be from 0 to pi                                                                                                  
# kappa_data: The density*opacity array (same size as theta)                                                                                            
# rho_data:   The density arrary (same size as theta)
# radius:     The radius                                                                                                                                   
# mode:       1: using sqrt( (kappa-kappaes)*kappaes) 2. using sqrt(kappa*kappaes) Here kappa can be 
#             either Rossland or Planck mean.
# hemi:       Which hemisphere we are integrating (upper or lower)   
######################################################################################################## 
#------------------------------------------------------------------------------------------------------#
# Aux Function: Eff_Kappa -- Calculate the effective kappa using sqrt( (kappa-kappaes)*kappaes)
#------------------------------------------------------------------------------------------------------#
def Eff_Kappa(kappaes, kappa_data, rho_data):                                                                                                           
    #print(kappa_data, rho_data, kappa_data/rho_data, kappaes)
    #print ( max(0.0, np.sqrt( kappaes*(kappa_data/rho_data - kappaes) )) )                                                                                  
    return np.sqrt( max(0.0, kappaes*(kappa_data/rho_data - kappaes)  )) 
#------------------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------------------#                                                                    
# Aux Function: Eff_Kappa2 -- Calculate the effective kappa using sqrt(kappa_p*kappaes)                                                                             
#------------------------------------------------------------------------------------------------------#                                                                    
def Eff_Kappa2(kappaes, kappa_data, rho_data):                                                                                                                               
    #print(kappa_data, rho_data, kappa_data/rho_data, kappaes)                                                                                                              
    #print ( max(0.0, np.sqrt( kappaes*(kappa_data/rho_data - kappaes) )) )                                                                                                 
    #print(kappaes, kappa_data, rho_data,  np.sqrt( max(0.0, kappaes*kappa_data  )))
    return np.sqrt( max(0.0, kappaes*kappa_data*rho_data  ))         
#------------------------------------------------------------------------------------------------------#


def Get_Photosphere2(theta_data, kappa_data, rho_data, radius, mode, hemi):
    data_len = len(theta_data) 

    if hemi == 'lower':                                                                                                                                 
        kappa_data = kappa_data[::-1]                                                                                                                           
        theta_data = theta_data[::-1]
        rho_data   = rho_data[::-1]
      
    iphot = 0
    tau = np.zeros(data_len//2)
    for it in range(1, data_len//2, 1):
        if mode == 1:
            tau[it] = tau[it-1] + 0.5*radius*(Eff_Kappa(kappaes_code, kappa_data[it-1], rho_data[it-1])
            + Eff_Kappa(kappaes_code, kappa_data[it], rho_data[it]) )*(theta_data[it] - theta_data[it-1]) 
        elif mode == 2:
            tau[it] = tau[it-1] + 0.5*radius*(Eff_Kappa2(kappaes_code, kappa_data[it-1], rho_data[it-1])                        
            + Eff_Kappa2(kappaes_code, kappa_data[it], rho_data[it]) )*(theta_data[it] - theta_data[it-1])
        if np.abs(tau[it]) >= 1:    # found photosphere, no need to go further 
            #print("Hit!!!", tau[it], tau[it-1], kappa_data[it-1]/rho_data[it-1], kappa_data[it]/rho_data[it])
            iphot = it if (tau[it]+tau[it-1]) < 2 else it-1
            break
        iphot = it

    #print("data_len:", data_len, iphot, hemi)
    if hemi == 'lower': 
        iphot = data_len - iphot -1
    return iphot 




################################################################################################################################
# Fucntion: Check_Load_Files  -- 1. For a range [start,end], check if the current iteration has been saved before
#                             -- 2. If yes, load the saved data and go to the first one that hasn't been saved
#                             -- 3. If all iters are saved, pass
#                             -- 4. If [start,end] doesn't overlap with saved files or no saved files at all, output new arrays
################################################################################################################################
# Parameters: 
################################################################################################################################
# filenames:  quant filenames in the checkpoint directory that we check for. They should contain start and end iters.  
# t_filenames_pre:  Prefix of time filenames in the checkpoint directory.
# file_exist:  A flag showing if we find a matching save files or not.
# start:   Assigned starting iteration to compare with saved files:
# end:   Assigned ending iteration to compare with saved files:
# output_all: Output everything? True (t, quant_data, file_exist, save_start, start) | False (quant_data only)
# read_time: Read in time file? Set it False if you only want to regenerate some quantity files
# succinct: If we print succinct ('True') info or detailed ('False') info.
################################################################################################################################

def Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, output_all, read_time, succinct):
    Print_text( "Function Check_Load_Files: assgined start iter:", start, " end iter:", end)
    save_file_list = glob.glob(filenames)

    if len(save_file_list) == 0:
        Print_subsubtitle("No previously saved files, start new arrays")                                                                                                                                            
        t = []                                                                                                                                                                                                      
        quant_data = []                                                                                                                                                                                             
        save_start = start      

    else:
        save_list = []
        if len(save_file_list) > 0:                                                                                          
            for file in save_file_list:                  
                save_start = int(file.split('_')[-3])                                                                         
                save_end = int(file.split('_')[-2])                                                                           
                save_list.append([save_start, save_end])

        if succinct == "False": Print_text( "Available files have start/end iterations:", save_list)


        max_idx = 0
        max_end = start
        for i in range(0, len(save_list)):            
            if save_list[i][0] >= save_start and save_list[i][0] <= save_end:                 
                file_exist = 1
                if save_list[i][1] >= max_end:
                    max_idx = i
                    max_end = save_list[i][1]
            else:
                Print_text("File", save_file_list[i], "not overlaps with assigned period [",start, end,"]") 

        if succinct == "False": Print_text("Found best file, ", max_idx, save_list[max_idx], save_file_list[max_idx])



        file = save_file_list[max_idx]
        save_start = save_list[max_idx][0]
        save_end = save_list[max_idx][1]     
 
        if succinct == "False": Print_text( "Function Check_Load_Files: saved start iter:", save_start, "saved end iter:", save_end)
       
        if file_exist == 1:
            if end <= save_end:                                                                                      
                Print_subsubtitle("The entire range",start, end, "is included in the file:", file)                                             
            else:                                                                                                    
                Print_text("Need to load and read new data. The file we load is:", file)                                                                   

            quant_data = np.load(file)                                                                               
            Print_text("Structure of the readin file:", np.shape(quant_data))

            if output_all:
                if read_time == "True":
                    t_file = t_filenames_pre + str(save_start)+"_"+str(save_end)+"_iter.npy"                                       
                    t = np.load(t_file)                                                                                      
                else:
                    t = []

                start = save_end + 1                 
                if start < end:
                   Print_subsubtitle("Move the start to preivous end. Now we load data from", start, "to", end)

            Print_text("Done reading in data, the dataset covers iteration:", save_start, " to:", save_end)                  
        else: #file_exist == 0  
            Print_subsubtitle("No saved file overlaps with assigned period, start new arrays")                                       
            t = []
            quant_data = []
            save_start = start            

    if output_all:
        if succinct == "False": Print_subsubtitle("Output all is True!!!") 
        return t, quant_data, file_exist, save_start, start
    else:
        if succinct == "False": Print_subsubtitle("Output all is False!")
        return quant_data


###############################################################################################################################
#Fucntion: Save_Files  -- Save data to files at given iterations. Delete previouly saved files. 
###############################################################################################################################
# Parameters:
###############################################################################################################################
# save_step:  Save to file every # iterations
# iter: Iteration number to save files (optional, set it to None when not used)
# save_start: Starting iteration of the data (optional, set it to None when not used)
# start: Starting iteration of the assigned period (needed to avoid index error)
# end: Ending iteration of the assigned period (needed to save the last iteration)
# quant_data: Data to save to files
# save_file_pre: Prefix of the saved filename
# succinct: If we print succinct ('True') info or detailed ('False') info
###############################################################################################################################

def Save_Files(save_step, iter, save_start, start, end, quant_data, save_file_pre, succinct):
    name_body = "_"
    for item in [str(save_start), str(iter)]:
        if item != "None":
            name_body += item
            name_body += "_"
 
    np.save(save_file_pre + name_body + "iter.npy", quant_data)  
    Print_subsubtitle("Function Save_Files:", save_file_pre + name_body + "iter.npy")
    if succinct == "False":
        Print_text("Now the file contains data with size:", np.shape(quant_data))
       
    if str(save_start) != "None" and str(iter) != "None":
        filename_pre = save_file_pre + "_" + str(save_start) + "_" + str(iter-save_step) + "_" + "iter.npy"
        if os.path.exists(filename_pre):
            if succinct == "False": 
                Print_subsubtitle("Deleting files", filename_pre)
            os.remove(filename_pre)



###############################################################################################################################
#Fucntion: Plotting_Mesh  -- Make the mesh plot of the data
###############################################################################################################################
# Parameters:
###############################################################################################################################
# x: x-axis quantity
# y: y-axis quantity
# q: Major quantity to plot
# x_fac: Factor to convert x-unit
# y_fac: Factor to convert y-unit
# q_fac: Factor to convert q-unit
# log_x: Set x in log scale? 
# log_y: Set y in log scale?
# log_q: Set q in log scale?
# color: Colormap used
# x_lim: x-limit (list)
# y_lim: y-limit (list) 
# v_max: q max-limit (can be None)
# v_min: q min-limit (can be None) 
# title: Figure title
# xlabel: Figure xlabel
# ylabel: Figure ylabel
# cbar_title: Colorbar title (optional)
# font_dict: the font used for title, label, etc.
###############################################################################################################################

def Plotting_Mesh(x,y,q,x_fac,y_fac,q_fac,log_x,log_y,log_q, color, x_lim, y_lim, v_max, v_min, title, xlabel, ylabel, cbar_title, font_dict):

    Print_subtitle("Function Plotting_Mesh: Starting making plot!")

    X, Y = np.meshgrid(x, y) 
    
    Print_subsubtitle(X, Y, np.shape(X), np.shape(Y))
    Print_subsubtitle(np.shape(q))   

    X *= x_fac
    Y *= y_fac
    if len(np.shape(q)) == 3:
       q = [ [qqq*q_fac for qqq in qq] for qq in q]
    elif len(np.shape(q)) == 2:
       q = [ qq*q_fac for qq in q] 
    
    if log_x:
        print("x in logscale!")
        X_plot = np.log10(X)
    else:
        X_plot = X

    if log_y:
        print("y in logscale!")                                                                                                        
        Y_plot = np.log10(Y) 
    else:                                                                                                           
        Y_plot = Y

    if log_q:   
        print("q in logscale!")                                                                                                     
        Q_plot = np.log10(q) 
    else:                                                                                                           
        Q_plot = q

    fig, ax1 = plt.subplots(1, 1)
 
    print(np.shape(X_plot), np.shape(Y_plot), np.shape(Q_plot))
    if v_min != None and v_max != None:
        im0 = ax1.pcolormesh(X_plot, Y_plot, Q_plot, cmap=color, shading='auto', vmax=v_max, vmin=v_min)    
    else:    
        im0 = ax1.pcolormesh(X_plot, Y_plot, Q_plot, cmap=color, shading='auto')

    cbar =fig.colorbar(im0, ax=ax1, fontdict = font_dict)
    if cbar_title != "None":
        cbar.ax.set_ylabel(cbar_title, rotation = 90)

    ax1.set_title(title, fontdict = font_dict)
    ax1.set_xlabel(xlabel, fontdict = font_dict)
    ax1.set_ylabel(ylabel, fontdict = font_dict)
   
    if x_lim != None:
        ax1.set_xlim(x_lim[0], x_lim[1])
    if y_lim != None:                                                                                                
        ax1.set_ylim(y_lim[0], y_lim[1])
   
    return ax1


###############################################################################################################################
#Fucntion: Plotting_Mesh_YX  -- Same as Plotting_Mesh but switch mesh order between X and Y
###############################################################################################################################
# Parameters:
###############################################################################################################################
# (same of Plotting_Mesh)
###############################################################################################################################

def Plotting_Mesh_YX(x,y,q,x_fac,y_fac,q_fac,log_x,log_y,log_q, color, x_lim, y_lim, v_max, v_min, title, xlabel, ylabel, cbar_title, font_dict):      

    Print_subtitle("Function Plotting_Mesh_YX: Starting making plot!")                                                                             
                                                                                                                                          
    Y, X = np.meshgrid(y, x)                                                                                                                                                                                       
                                                                                                                                                                                                                    
    Print_subsubtitle("Data Structure of X, Y data: ", np.shape(X), np.shape(Y))                                                                                                                                                               
    Print_subsubtitle(np.shape(q)) 
                                                                                                                                                                                                                    
    X *= x_fac                                                                                                                                                                                                      
    Y *= y_fac                                                                                                                                                                                                      
    q = [ [qqq*q_fac for qqq in qq] for qq in q]                                                                                                                                                                    
                                                                                                                                                                                                                    
                                                                                                                                                                                                                    
    if log_x:                                                                                                                                                                                                       
        print("x in logscale!")                                                                                                                                                                                     
        X_plot = np.log10(X)                                                                                                                                                                                        
    else:                                                                                                                                                                                                           
        X_plot = X                                                                                                                                                                                                  
                                                                                                                                                                                                                    
    if log_y:                                                                                                                                                                                                       
        print("y in logscale!")                                                                                                                                                                                     
        Y_plot = np.log10(Y)                                                                                                                                                                                        
    else:                                                                                                                                                                                                           
        Y_plot = Y                                                                                                                                                                                                  
                                                                                                                                                                                                                    
    if log_q:                                                                                                                                                                                                       
        print("q in logscale!")                                                                                                                                                                                     
        Q_plot = np.log10(q)                                                                                                                                                                                        
    else:                                                                                                                                                                                                           
        Q_plot = q                                                                                                                                                                                                  
                                                                                                                                                                                                                    
    fig, ax1 = plt.subplots(1, 1)                                                                                                                                                                                   
                                                                                                                                                                                                                    
    print(np.shape(X_plot), np.shape(Y_plot), np.shape(Q_plot))                                                                                                                                                     
    if v_min != None and v_max != None:                                                                                                                                                                             
        im0 = ax1.pcolormesh(X_plot, Y_plot, Q_plot, cmap=color, shading='auto', vmax=v_max, vmin=v_min)                                                                                                            
    else:                                                                                                                                                                                                           
        im0 = ax1.pcolormesh(X_plot, Y_plot, Q_plot, cmap=color, shading='auto')                                                                                                                                    
                                                                                                                                                                                                                    
    cbar =fig.colorbar(im0, ax=ax1)                                                                                         
    if cbar_title != "None":                                                                                                
        cbar.ax.set_ylabel(cbar_title, rotation = 90, fontdict = font_dict) 
                                                                                                                                                                                                                    
    ax1.set_title(title, fontdict = font_dict)                                                                                                                                                                      
    ax1.set_xlabel(xlabel, fontdict = font_dict)                                                                                                                                                                    
    ax1.set_ylabel(ylabel, fontdict = font_dict)                                                                                                                                                                    
                                                                                                                                                                                                                    
    if x_lim != None:                                                                                                                                                                                                
        ax1.set_xlim(x_lim[0], x_lim[1])
    if y_lim != None:       
        ax1.set_ylim(y_lim[0], y_lim[1])
                                                                                                                                                                                                                    
    return ax1                                                                                                                                                                                                      

###############################################################################################################################
#Fucntion: Plotting_Mesh_MeshXY  -- Same as Plotting_Mesh but with X/Y already meshed                                                                                                                        
###############################################################################################################################
# Parameters:                                                                                                                                                                                                       
###############################################################################################################################
# (same of Plotting_Mesh)                                                                                                                                        
###############################################################################################################################
                                                                                                                                                                                                                    
def Plotting_Mesh_MeshXY(x,y,q,x_fac,y_fac,q_fac,log_x,log_y,log_q, color, x_lim, y_lim, v_max, v_min, title, xlabel, ylabel, cbar_title, font_dict):     
                                                                                                                                                                                                                    
    Print_subtitle("Function Plotting_Mesh_MeshXY: Starting making plot!")

    X = [ [xxx*x_fac for xxx in xx] for xx in x]                                                                                                                                                                                                      
    Y = [ [yyy*y_fac for yyy in yy] for yy in y]                                                                                                                                                                                                      
    q = [ [qqq*q_fac for qqq in qq] for qq in q]                                                                                                                                                                    
                                                                                                                                                                                                                    
    Print_subsubtitle("Data Structure of X, Y data: ", np.shape(X), np.shape(Y))                                                                                                                                                               
    Print_subsubtitle(np.shape(q)) 

                                                                                                                                                                                                                    
    if log_x:                                                                                                                                                                                                       
        print("x in logscale!")                                                                                                                                                                                     
        X_plot = np.log10(X)                                                                                                                                                                                        
    else:                                                                                                                                                                                                           
        X_plot = X                                                                                                                                                                                                  
                                                                                                                                                                                                                    
    if log_y:                                                                                                                                                                                                       
        print("y in logscale!")                                                                                                                                                                                     
        Y_plot = np.log10(Y)                                                                                                                                                                                        
    else:                                                                                                                                                                                                           
        Y_plot = Y                                                                                                                                                                                                  
                                                                                                                                                                                                                    
    if log_q:                                                                                                                                                                                                       
        print("q in logscale!")                                                                                                                                                                                     
        Q_plot = np.log10(q)                                                                                                                                                                                        
    else:                                                                                                                                                                                                           
        Q_plot = q                                                                                                                                                                                                  
                                                                                                                                                                                                                    
    fig, ax1 = plt.subplots(1, 1)                                                                                                                                                                                   
                                                                                                                                                                                                                    
    print(np.shape(X_plot), np.shape(Y_plot), np.shape(Q_plot))                                                                                                                                                     
    if v_min != None and v_max != None:                                                                                                                                                                             
        im0 = ax1.pcolormesh(X_plot, Y_plot, Q_plot, cmap=color, shading='auto', vmax=v_max, vmin=v_min)                                                                                                            
    else:                                                                                                                                                                                                           
        im0 = ax1.pcolormesh(X_plot, Y_plot, Q_plot, cmap=color, shading='auto')                                                                                                                                    
                                                                                                                                                                                                                    

    cbar =fig.colorbar(im0, ax=ax1)                                                                                         
    if cbar_title != "None":                                                                                                
        cbar.ax.set_ylabel(cbar_title, rotation = 90, fontdict = font_dict) 
                                                                                                                                                                                                        
    ax1.set_title(title, fontdict = font_dict)                                                                                                                                                                      
    ax1.set_xlabel(xlabel, fontdict = font_dict)                                                                                                                                                                    
    ax1.set_ylabel(ylabel, fontdict = font_dict)                                                                                                                                                                    
                                                                                                                                                                                                                    
    if x_lim != None:                                                                                                                                                                                                
        ax1.set_xlim(x_lim[0], x_lim[1])
    if y_lim != None:       
        ax1.set_ylim(y_lim[0], y_lim[1])
                                                                                                                                                                                                                    
    return ax1                                                                                                                                                                                                      
                                                                                                                                                                                                                    
###############################################################################################################################
#Fucntion: Plotting_1D  -- Make 1D plot on top of the mesh plot
###############################################################################################################################
# Parameters:  
###############################################################################################################################
# ax: ax object from the mesh plot
# x: x-axis quantity
# x_fac: Factor to convert x-unit
# log_x: Set x in log scale? 
# quant_list: List of quantites to plot
# q_fac: Factor to convert q-unit   
# log_q: Set q in log scale?
# plt_color: Color
# xlabel: x-axis label
# ylabel: y-axis label
# label_q: plotting legend
# font_d:  font dictionary
###############################################################################################################################
def Plotting_1D(ax, x, x_fac, log_x, quant_list, q_fac, log_q, plt_color, xlabel, ylabel, label_q, font_d):
    x = [xx * x_fac for xx in x]
    if log_x:
        x = np.log10(x)

    if len(np.shape(quant_list)) == 2:
        for item in quant_list:
            item = [i * q_fac for i in item]
            if log_q:
                item = np.log10(item)  
            ax.plot(x, item, color = plt_color, label = label_q)
            if label_q != None:
                ax.legend()
            if xlabel != None:
                ax.set_xlabel(xlabel, fontdict = font_d)
            if ylabel != None:
                ax.set_ylabel(ylabel, fontdict = font_d)  
    elif len(np.shape(quant_list)) == 1:  
        item = [ i * q_fac for i in quant_list] 
        if log_q:
            item = np.log10(item)
        ax.plot(x, item, color = plt_color, label = label_q)
        if label_q != None:  
            ax.legend()
        if xlabel != None:
            ax.set_xlabel(xlabel, fontdict = font_d)
        if ylabel != None:
            ax.set_ylabel(ylabel, fontdict = font_d)
    else:
        print("Wrong y value structure!!!")
        os.exit()





