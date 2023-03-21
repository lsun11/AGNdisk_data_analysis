# Welcome to use the AGNdisk Data Analysis code for <span style="font-family: 'Helvetica', sans-serif;">Atheta++</span>
![alt text](https://github.com/lsun11/AGNdisk_data_analysis/blob/main/Data_Analysis_Code/readme_pics/code_screen_shot.png)
(↑ The screen shot of the starting interaface of the code!)

## Introduction

This is a codebase mainly consists of Python code controlled by bash commands scripts that can automatically create beautiful and professional plots by loading/reading the output `.npz` files from the <span style="font-family: 'Helvetica', sans-serif;"> Athena++</span> astrophysical magnetohydrodynamics (MHD) code in C++

The data files code can read are in `.npz` format. Each file contains arrays(1D or 2D) of physical quantities. For exmample:

     >>> data = np.load('../DATAWedge8/hist_02000.npz')
     >>> data.files
     ['time', 'radius', 'rloc', 'theta', 'x1f', 'x2f', 'rho', 'Er', 'sigma', 'sigma_p', 'Fr1', 'Fr2', 'B1', 'B2', 'B3', 'B1sq', 'B2sq', 'B3sq', 'PB1', 'PB2', 'PB3', 'pg', 'rhovr', 'MRIlambda', 'Ek1', 'Ek2', 'Ek3', 'BrBp', 'BtBp', 'rhovrvphi', 'rhovtvphi', 'rhovp', 'rhovt']

The output array aboves contains headers of the physical quantities. The files are named based on iteration (time) (e.g. `hist_02000.npz`) and are stored in directory `/DATAWedge8`, which is named based on the simulation cases (Wedge8, Wedge9B, Wedge10). When launching the code, the data will be mounted from the remote sever (KITP) to the local directory based on the cases to analyze. The mounting process requires a password. 


## How to use:  
### 1. Before using the code, check the required packages and the remote connection.

You need the following 3 packages to run the code:  
&nbsp; &nbsp; &nbsp; 1). figlet https://linuxhint.com/figlet-command-linux

<font size="2">(note: some version of the figlet may not have the font in the screenshot above, you can simply switch it to another font in `update.sh` or just leave it blank)</font>

&nbsp; &nbsp; &nbsp; 2). colorama https://pypi.org/project/colorama

&nbsp; &nbsp; &nbsp; 3). sshfs https://www.digitalocean.com/community/tutorials/how-to-use-sshfs-to-mount-remote-file-systems-over-ssh  

You also need to check if your working environment is connected to the remote data base (e.g. KITP server):  
`ssh #username@128.111.9.179`


### 2. Download the code in the location directory and run a quick test.  
You can run a (relatively) quick test by hitting `./update_test.sh` in the `/Command` directory. This test plots a very short period of all types of data and it takes about 40 mintues. For full run some data requries reading huge amount of data and could take hours or even a day to finish. 

&nbsp; &nbsp;  &nbsp;  &nbsp; **Note**: It's recommended to rename the plotting and checkpoint directories in parameter.py before running the test:

 i.e. change: 
 
       Plot_dir             = Data_dir + 'PLOTS/'
       Checkpoint_dir       = Data_dir + 'Checkpoint/'
       
  to something like:
      
       Plot_dir             = Data_dir + 'PLOTS_Test/'
       Checkpoint_dir       = Data_dir + 'Checkpoint_Test/'
so the test plots and checkpoints will be saved in separated places. For production run you can change it back or rename to whatever you like.   

If the test runs successfully, you should see plots are saved in `Plot_dir` you set (see the following picture as an example) and checkpoints files in `Checkpoint_dir`:

![alt text](https://github.com/lsun11/AGNdisk_data_analysis/blob/main/Data_Analysis_Code/readme_pics/test.png)




### 3. To run the code in full length, you need to modify the script `input.sh` in the `/Command` directory.

&nbsp; &nbsp; &nbsp; **1). The first** section in `input.sh` has variables that needed to be modified for each run.

  `case` : simulation case you want to plot. If is this set as a list (e.g.  (Wedge8 Wedge9B Wedge10) ), each case in the list will be plotted.
  `radii`: many plots requires data at a given radius, set the radius here. This can also in form of a list (e.g. (120 140 160) ) for multiple runs.
  `start` and `end`: starting and ending iteration you want to analyze. The code will automatically compare your input start and end with the available
  iterations. If one or both of them are out of range, the code will auto set to the maximum possible value or instruct to reset. 


&nbsp; &nbsp; &nbsp; **2). The second** section is the main body of the command. It controls wether or not to make each types of plot and some specific parameters of making the plot. Current the code supports the following types of plots:

   #### &nbsp; &nbsp; &nbsp;Make poloidal distribution plot of azimuthally-averaged quantity

   #### &nbsp; &nbsp; &nbsp;Make quantities vs times vs theta (with photosphere) (quantities e.g. density, $T_{rad}$, $T_{gas}$, ($S_h^{r\phi}-S_m^{r\phi}$), ($S_h^{\theta\phi}-S_m^{\theta\phi}$) )
   
   
   #### &nbsp; &nbsp; &nbsp;Make time-average quantities vs theta or time   (e.g. $P_B$, $E^\theta_k$, $P_{gas}$, $P_{rad}$)
   
   #### &nbsp; &nbsp; &nbsp;Make theta-average quantities vs r   (Most likely $\dot{M}$) 
   
   #### &nbsp; &nbsp; &nbsp;Make quantity vs time vs radius plots at a given theta 
   
   #### &nbsp; &nbsp; &nbsp;Make temperature variation plots at given locations
   
For each case, the first parameter is the switch of making this type of the plot or not. Set Y/N for make plot/don't make plot. Note that in case **Make quantities vs times vs theta (with photosphere)**, each quantity has its own swith. 

The other parameters for each varies but there are several common ones with the same suffix:

`_plot_phot`: compute and plot photoshpere and plot on top? Y/N  
`_phot_mode`: method to the compute the photosphere.   
        
&nbsp; &nbsp; &nbsp;-1: use the absolute value.  
        
&nbsp; &nbsp; &nbsp;  1: use $\sqrt{ (\kappa-\kappa_{es})*\kappa_{es}}$
        
&nbsp; &nbsp; &nbsp;  2: use $\sqrt{\kappa * \kappa_{es}}$  
&nbsp; &nbsp; &nbsp; (Here $\kappa$ can be either Rossland or Planck.)  
 
`_save_t`   : save/load the time data or not? Y/N  
`_succ`     : do we output succinctly? Set it to N for debugging (or just to want to know the process better) Y/N  
`_rad_cut`  : upper boundary of the radius to plot, this is in unit of index.   


&nbsp; &nbsp; &nbsp; **3). The third** section has directories parameters and username for monting the data. These should be set correctly when use the code for the fist time. The default setting should give a correct environment to save/load data and make plots. 

&nbsp; &nbsp; &nbsp; **4). The fourth** section contains parameters that controls the appearance of the interface output. You are welcomed to play with it to make it look better! 

### 4. The analyzing/plotting procedure.

Once the `input.sh` file is set properly, we can start launching the code by typing `./update.sh`.

After the starting screen (shown above), the first stage is **initialization**:
![alt text](https://github.com/lsun11/AGNdisk_data_analysis/blob/main/Data_Analysis_Code/readme_pics/1.png)

Here the code checks if all the directories are set properly. If the directories don't exsit, the code creates them. If they exisit, the code shows the contents inside.


The second stage is **data mounting**:
![alt text](https://github.com/lsun11/AGNdisk_data_analysis/blob/main/Data_Analysis_Code/readme_pics/2.png)

Here the code would link the remote database to local directories if the data is not already linked. You need input the password of your remote access point here:
![alt text](https://github.com/lsun11/AGNdisk_data_analysis/blob/main/Data_Analysis_Code/readme_pics/2-2.png)


The third stage is **plotting**:
![alt text](https://github.com/lsun11/AGNdisk_data_analysis/blob/main/Data_Analysis_Code/readme_pics/3.png)

Now the code starts read/load data and make plots based on input.sh. More details will output if the succtint parameter is set to "False", otherwise the output are mostly in forms of progress bars.  


Lastly, the entire process is completed if the following appears:
![alt text](https://github.com/lsun11/AGNdisk_data_analysis/blob/main/Data_Analysis_Code/readme_pics/4.png)



