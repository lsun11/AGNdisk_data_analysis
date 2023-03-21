# Welcome to use the AGNdisk Data Analysis code for Atheta++
![alt text](https://github.com/lsun11/AGNdisk_data_analysis/blob/main/Data_Analysis_Code/code_screen_shot.png)

This is a codebase mainly consists of Python code controlled by bash commands scripts that can automatically create beautiful and professional plots by loading/reading the output .npz files from the Athena++ astrophysical magnetohydrodynamics (MHD) code in C++

How to use:
1. Before using the code, you need check if the following packages are installed:

&nbsp; &nbsp; &nbsp; 1). figlet https://linuxhint.com/figlet-command-linux/

&nbsp; &nbsp; &nbsp; 2). colorama https://pypi.org/project/colorama/

2. Download the code in the location directory. First you can run a (relatively) quick tests by hitting ./update_test.sh in the /Command directory. This tests plots a very short period of all types of data and it takes about 20 mintues. For full run some data requries reading huge amount of data and could take hours or even days to finish. 

&nbsp; &nbsp;  &nbsp;  &nbsp; **Note**: It's recommended to rename the plotting and checkpoint directories in parameter.py before running the test:

 i.e. Change: 
 
       Plot_dir             = Data_dir + 'PLOTS/'
       Checkpoint_dir       = Data_dir + 'Checkpoint/'
       
 To something like:
      
       Plot_dir             = Data_dir + 'PLOTS_Test/'
       Checkpoint_dir       = Data_dir + 'Checkpoint_Test/'
so the test plots and checkpoints will be saved in separated places. For production run you can change it back or rename to whatever you like. 
     
