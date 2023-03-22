############################### Stuff to change for each run ###############################
### cases, times, and radii ###
#cases=(Wedge8 Wedge9B Wedge10)
cases=(Wedge8)
radii=(120 140)
start=2000
end=2030

### type plots  you want to make ###
#plot_Mdot=Y

############################### Different types of plot to make ###############################
### Make poloidal distribution plot of azimuthally-averaged quantity
plot_sin_cos_quant=Y

plot_sin_cos_q1=rho
plot_sin_cos_q2=sigma_p
plot_sin_cos_iter=(2000 2400)
plot_sin_cos_plot_phot=Y
plot_sin_cos_phot_mode=-1
plot_sin_cos_succ=True


### Make density vs times vs theta (with photosphere)
plot_t_th=Y
plot_t_th_q=(rho B2 B3 Temp_r Temp_g Ang_rp Ang_tp PB2 Er)

plot_t_th_kappa=Y
plot_t_th_kappa_p=Y 

plot_t_th_sigma=sigma  #(sigma/sigma_p)
plot_t_th_phot=Y
plot_t_th_phot_mode=-1
plot_t_th_phot_save_t=True
plot_t_th_phot_succ=True




### Make time-average quantities vs theta or time  
plot_1D_tave=Y
plot_1D_tave_vst=Y

plot_1D_tave_q1=PB_total
plot_1D_tave_q2=pg
plot_1D_tave_q3=Ek2
plot_1D_tave_q4=Er
plot_1D_tave_save_t=True
plot_1D_tave_succ=True


### Make theta-average quantities vs r
plot_1D_th_ave=Y

plot_1D_th_ave_q1=rhovr
plot_1D_th_ave_q2=None
plot_1D_th_ave_q3=None
plot_1D_th_ave_q4=None
plot_1D_th_ave_Mdot=True
plot_1D_th_ave_save_t=True
plot_1D_th_ave_succ=True


### Make quantity vs time vs radius plots at a given theta
plot_t_r_quant=Y

plot_t_r_q1=sigma_p
plot_t_r_q2=None
plot_t_r_theta=90.0
plot_t_r_succ=True


### Make temperature variation plots at given locations
plot_T_var=Y

plot_T_var_deg=photosphere
plot_T_var_q1=Er
plot_T_var_q2=rho
plot_T_var_q3=sigma
plot_T_var_ph_mode=-1
plot_T_var_rad_cut=200
plot_T_var_succ=True

########################### Stuff to modify when initialzie the code #########################


### path to store the code ###
#home_path="/home/lunan/mnt/AGNWedge8_product2/"
home_path=$(dirname $(dirname $(pwd)))"/"
code_path=$home_path"Data_Analysis_Code/"
plot_path=$home_path"PLOTS_TEST/"
checkpoint_path=$home_path"Checkpoint_TEST/" 
username=luan




############################## Aux stuff (no need to change) ###############################
### 1. banner and boarders ###
# Define the banner text
banner_text="AGN ANALYSIS/PLOTTING CODE"
author_text="Maintained by: Lunan Sun"
# Get the width of the terminal window   
terminal_width=$(tput cols)
# Calculate the center position for the banner text
center_position=$((($terminal_width-${#banner_text})/2))
# Set the font size  
font_size=7

# set the length and char of the boarders
len=$terminal_width 
len2=$center_position
len3=$((($terminal_width-13)/2))
len4=$((($terminal_width-12)/2))
ch='#'
ch2='='
ch3='-'
ch4='|'
ch5='~'
thick_b=$(printf '%*s' "$len" | tr ' ' "$ch")
thick_b_short=$(printf '%*s' "$len2" | tr ' ' "$ch")
mid_b=$(printf '%*s' "$len" | tr ' ' "$ch2")
mid_b_short=$(printf '%*s' "$len3" | tr ' ' "$ch2")
thin_b=$(printf '%*s' "$len" | tr ' ' "$ch3")
thin_b_short=$(printf '%*s' "$len4" | tr ' ' "$ch3")
vert_b=$(printf '%*s' "$len" | tr ' ' "$ch4")
vert_b_short=$(printf '%*s' "$len3" | tr ' ' "$ch4")
wave_b=$(printf '%*s' "$len" | tr ' ' "$ch5")
wave_b_short=$(printf '%*s' "$len4" | tr ' ' "$ch5")

### 2. Date and time ###
current_time=$(date +"%Y-%m-%d %H:%M:%S")
