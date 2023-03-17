############################### Stuff to change for each run ###############################
### cases and radii ###
#cases=(Wedge8 Wedge9B Wedge10)
cases=(Wedge8)
radii=(120 140)

### type plots  you want to make ###
#plot_Mdot=Y

############################### Different types of plot to make ###############################
### Make poloidal distribution plot of azimuthally-averaged quantity
plot_sin_cos_quant=N
plot_sin_cos_q1=rho
plot_sin_cos_q2=sigma_p
plot_sin_cos_iter=(2000 2400)
plot_sin_cos_plot_phot=Y
plot_sin_cos_phot_mode=-1
plot_sin_cos_succ=False


### Make density vs times vs theta (with photosphere)
plot_t_th_dens=Y
plot_t_th_sigma=sigma  #(sigma/sigma_p)
plot_t_th_phot=Y
plot_t_th_phot_mode=-1
plot_t_th_phot_save_t=True
plot_t_th_phot_succ=True

########################### Stuff to modify when initialzie the code #########################


### path to store the code ###
home_path="/home/lunan/mnt/AGNWedge8_product2/"
code_path=$home_path"Data_Analysis_Code/"
plot_path=$home_path"PLOTS/"
checkpoint_path=$home_path"Checkpoint/" 
username=luan




############################## Aux stuff (no need to change) ###############################
### 1. banner and boarders ###
# Define the banner text
banner_text="AGN ANALYSIS/PLOTTING CODE"
author_text="Author: Lunan Sun"
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
