#THIS SCRIPT PLOTS THE TEMPERATURE VARIATIONS VS TIME AND RADIUS
#Location can be chosen at: midplane, theta = 80deg, theta = 70deg, and photosphere
#Example command: python plotting_T_var.py Wedge8_2 Er 3000 3050 photosphere rho kappa False

print("========================== Starting plotting_T_var.py =============================")
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
verbose = argv[8]
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
print(filename)                                                                                                      
data = np.load(dir+'/'+filename)                         

r = Get_All_1D('radius', data, 1100, dir)   
if deg == 'photosphere':
    th = Get_All_1D('theta', data, -1, dir)
########################################################################################################  



if str(case)[-2:] == '_2':                                                                                                                                                                                                        
    print('case: _2')                                                                                                                                                                                                             
    idx = {"surface_density":6, 'rho':7, 'Er':8, 'kappa':9, 'Fr1':10, 'Fr2':11, 'B1_r1':12, 'B1_r2':13, 'B1_r3':14, 'B1_r4':15, 'B2_r1':16, 'B2_r2':17, 'B2_r3':18, 'B2_r4':19, 'B3_r1':20, 'B3_r2':21, 'B3_r3':22, 'B3_r4':23, 'PB_r2':25, 'PB_r3':26, 'PB_r4':27, 'pg_r1':28, 'pg_r2':29, 'pg_r3':30, 'pg_r4':31, 'rhovr':32, 'lambda_r1':33, 'lambda_r2':34, 'lambda_r3':35, 'lambda_r4':36}
else:                                                                                                                                                                                                       
    print('case: no _2')                                                                                                                                                                                                          
    idx = {"surface_density":6, 'rho1':7, 'rho2':8, 'rho3':9, 'rho4':10, 'Er1':11, 'Er2':12, 'Er3':13, 'Er4':14, 'kappa1':15, 'kappa2':16, 'kappa3':17, 'kappa4':18} 


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
    
    t, quant_data, file_exist, save_start, start = Check_Load_Files(filenames, t_filenames_pre, file_exist, start, end, True, True, verbose)
######################################################################################################## 
#quant_data = quant_data[:-2]
print(np.shape(t), np.shape(quant_data), file_exist)





########################################################################################################
# Iterate from start to end to compute and save time and quant
######################################################################################################## 
for iter in range(start, end):                                                                                                                                                                                                
        filename = "hist_"+str(iter).zfill(5)+".npz"                                                                                                                                                                              
        print(filename)     
        data = np.load(dir+'/'+filename)      

        if file_exist == 1:
            t = np.append(t, Get_Time(data, dir))
        else:
            t.append(Get_Time(data, dir))



        if deg != "photosphere":
            print("We are plotting quantities at "+str(deg)+"!!!") 
            quant_data.append(data[quant][deg_idx][0:1100]) #midplane

        else:   #integrate rho*kappa to get photosphere
            print("We are plotting quantities at the photosphere!!!")

            r_data    = r[0:1100]                               # 0:1100 --> we only plot radius < 500 rg             
            kappa_zip = list(zip(*data[quant3]))
            kappa_data = kappa_zip[0:1100]      
            th_data   = th        
          
            print (len(kappa_data), len(kappa_data[0]))
            # initialize an array for the indices of photosphere theta
            iphot = []
            hemi = 'upper'

            for it_r in range(len(r_data)):
                iphot.append(Get_Photosphere(th_data, kappa_data[it_r], r_data[it_r], hemi))

            
            # This line below is SUPER SLOW! 
            if file_exist == 0:
                quant_data.append([data[quant][iphot[i]][i] for i in range(len(iphot))])
            else:
                quant_data = np.vstack([quant_data, [data[quant][iphot[i]][i] for i in range(len(iphot))]])
                print(np.shape(quant_data))

######################################################################################################## 
# Save files since appending above is too slow!!!
######################################################################################################## 
            save_step = 30
           
            if iter > start and (iter%save_step == 0 or iter == end-1):             
                save_quant_file_pre = checkpoint_path + "T_var_Er"
                Save_Files(save_step, iter, save_start, start, end, quant_data, save_quant_file_pre)

                save_time_file_pre = checkpoint_path + "T_var_time"
                Save_Files(save_step, iter, save_start, start, end, t, save_time_file_pre) 

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
print(np.shape(T_data))
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
xlabel = 'Time (s)'
ylabel = '$\log_{10}(r/r_g)$'

time_start = int(t[0]*time_to_sec)
savename = save_path+'photosphere_rad_temperature_var_'+str(case)+'_'+str(deg)+'_t_'+str(time_start)+'.png'
font_dict = font

Plotting_Mesh(x,y,q,x_fac,y_fac,q_fac,log_x,log_y,log_q, color, xlim, ylim, v_max, v_min, title, xlabel, ylabel, font_dict)
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

