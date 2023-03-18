# load the input file
source input.sh

# Print the banner (require to install figlet)
echo $(tput setaf 24)$thick_b 
echo $(tput setaf 24)$thick_b
echo $(tput setaf 24)$thick_b 
#### Somehow the figlet in eddington doesn't have font univers so I use block here...
#echo "$(tput setaf 53)$(tput bold)$(printf '%*s' $center_position)$(figlet -f univers -w $((terminal_width-$font_size)) "${banner_text}")$(printf '%*s' $center_position)$(tput sgr0)"
echo "$(tput setaf 160)$(tput bold)$(printf '%*s' $center_position)$(figlet -f block -w $((terminal_width-$font_size)) "${banner_text}")$(printf '%*s' $center_position)$(tput sgr0)" 
echo "$(tput setaf 17)$(tput bold)$(printf '%*s' $((terminal_width-${#author_text}-17)))${author_text}$(printf '%*s' ${#author_text} )$(tput sgr0)"
echo "$(tput setaf 17)$(tput bold)$(printf '%*s' $((terminal_width-${#current_time}-15)))"Current Time: "${current_time}$(printf '%*s' ${#current_time} )$(tput sgr0)"
echo $(tput setaf 24)$thick_b
echo $(tput setaf 24)$thick_b
echo $(tput setaf 24)$thick_b  
echo ""
echo "" 
sleep 1

echo $(tput setaf 6)$thick_b                                                                                                                                                                                                  
echo "$(figlet -f standard -w $((terminal_width-$font_size)) "1.  Initialize!!!")"          
echo $(tput setaf 6)$thick_b                                                                    
tput setaf 0 
./initialize.sh
sleep 1

echo ""
echo ""
echo $(tput setaf 6)$thick_b                                                                                                                                                                                                                   
echo "$(figlet -f standard -w $((terminal_width-$font_size)) "2.  Mount Data!!!")" 
echo $(tput setaf 6)$thick_b 
tput setaf 0
./mount.sh

echo ""
echo ""
echo $(tput setaf 6)$thick_b 
echo "$(figlet -f standard -w $((terminal_width-$font_size)) "3.  Start Making Plots!!!")"
echo $(tput setaf 6)$thick_b                                                                                                                                                                                                  
sleep 1

#####################################################
######### Start looping through cases ###############
#####################################################
tput setaf 0
for case in ${cases[@]} ; do
    echo ""

    tput setaf 24
    if [[ $case == "Wedge8" ]]; then  
        echo $mid_b
	echo $vert_b_short" Case Wedge8 "$vert_b_short
        echo $mid_b
	start_all=2000
	end_all=3658
    elif [[ $case == "Wedge8_2" ]]; then
      	echo $mid_b
        echo $vert_b_short" Case Wedge8_2 "$vert_b_short
	echo $mid_b
	start_all=2000
	end_all=3658
    elif [[ $case == "Wedge9B" ]]; then
        echo $mid_b
        echo $vert_b_short" Case Wedge9B "$vert_b_short
        echo $mid_b
	start_all=0
	end_all=1157
    elif [[ $case == "Wedge9B_2" ]]; then
        echo $mid_b
        echo $vert_b_short" Case Wedge9B_2 "$vert_b_short
        echo $mid_b
	start_all=0
	end_all=1213
    else
	echo "Error"
	exit 1;
    fi
    tput setaf 0

    ### determine which start and end iteration we use ###
    if [[ $start == "None" ]]; then 
        echo "We use earliest possible starting iteration: $start_all"
        start=$start_all
    elif [[ $start < $start_all ]]; then
        echo "Input starting iteration is smaller than possible range, setting start to earliest possible starting iteration: $start_all"
        start=$start_all
    elif [[ $start > $end_all ]]; then 
        echo "Error. Input starting iteration is larger than maximum possible iteration, please reset!!!"
        exit 1;        
    fi   
 
    if [[ $end == "None" ]]; then                                                                                                                                        
        echo "We use last possible ending iteration: $end_all"                                                                                                     
        end=$end_all   
    elif [[ $end > $end_all ]]; then
        echo "Input ending iteration is larger than possible range, setting end to last possible ending iteration: $end_all" 
        end=end_all 
    elif [[ $end < $start_all ]]; then
        echo "Error. Input ending iteration is smaller than minimum possible iteration, please reset!!!"                                                                  
        exit 1;  
    fi
 
    tput setaf 24
    echo $thin_b
    echo "$(printf '%*s' $(($center_position-20))) The starting and ending iterations we use are: $start, $end"
    echo $thin_b 
    tput setaf 0 

    #############################################
    ### Start plotting for each case ###
    ############################################# 
    if [[ $plot_sin_cos_quant =~ ^[Yy]$ ]]; 
    then
        tput setaf 31
        echo $thin_b
        echo "$(printf '%*s' $(($center_position-22)) ) Plotting poloidal distribution plot of azimuthally-averaged quantity! "
        echo $thin_b
        tput setaf 0
        for iter in ${plot_sin_cos_iter[@]}; do
            python $code_path"plotting_costh_vs_sinth_vs_quant.py" $case $plot_sin_cos_q1 $plot_sin_cos_q2 $iter $plot_sin_cos_plot_phot $plot_sin_cos_phot_mode $plot_sin_cos_succ;
        done
    fi

                                                                                                                                                                                                                                                    
    if [[ $plot_1D_th_ave =~ ^[Yy]$ ]];                                                                                                                                                                                                         
    then                                                                                                                                                                                                                                        
       tput setaf 31                                                                                                                                                                                                                            
       echo $thin_b    
       echo "$(printf '%*s' $(($center_position-11)) ) Plotting theta-average quantities vs radius    "                                                                                                                                         
       echo $thin_b                                                                                                                                                                                                                             
       tput setaf 0                                                                                                                                 
       python $code_path"plotting_1D_tave_vs_r.py" $case $start $end $plot_1D_th_ave_q1 $plot_1D_th_ave_q2 $plot_1D_th_ave_q3 $plot_1D_th_ave_q4 $plot_1D_th_ave_Mdot $plot_1D_th_ave_save_t $plot_1D_th_ave_succ
    fi 

    
    if [[ $plot_t_r_quant =~ ^[Yy]$ ]];                                                                                      
    then     
       tput setaf 31                                                                                                         
       echo $thin_b                                                                                                          
       echo "$(printf '%*s' $(($center_position-15)) ) Plotting quantity vs time vs radius plots at a given theta     " 
       echo $thin_b                                                                                                          
       tput setaf 0                                                                                                          
       python $code_path"plotting_t_vs_r_vs_quant.py" $case $plot_t_r_q1 $plot_t_r_q2 $plot_t_r_theta $start $end $plot_t_r_succ
    fi

    

    for radius in ${radii[@]}; do	
        tput setaf 24
        echo $thin_b
	echo $wave_b_short" radius=$radius "$wave_b_short
        echo $thin_b
        tput setaf 0                 

        if [[ $plot_t_th_dens =~ ^[Yy]$ ]];
        then
           tput setaf 31
           echo $thin_b 
           echo "$(printf '%*s' $(($center_position-13)) ) Plotting density vs time vs theta (with photosphere)   "
           echo $thin_b 
           tput setaf 0
           python $code_path"plotting_t_vs_th_vs_quant.py" $case rho $plot_t_th_sigma None $radius Y $start $end $plot_t_th_phot_mode $plot_t_th_phot_save_t $plot_t_th_phot_succ
        fi

        if [[ $plot_1D_tave =~ ^[Yy]$ ]];                                                                                                                                                                         
        then                                                                                                                                                                                                        
           tput setaf 31                                                                                                                                                                                            
           echo $thin_b                                                                                                                                                                                             
           echo "$(printf '%*s' $(($center_position-10)) ) Plotting time-average quantities vs theta   " 
           echo $thin_b                                                                                                                                                                                                        
           tput setaf 0       
           python $code_path"plotting_1D_tave_vs_theta.py" $case $start $end $radius $plot_1D_tave_q1 $plot_1D_tave_q2 $plot_1D_tave_q3 $plot_1D_tave_q4 vs_theta $plot_1D_tave_save_t $plot_1D_tave_succ
        fi
   
        if [[ $plot_1D_tave_vst =~ ^[Yy]$ ]];            
        then                                                                                                          
           tput setaf 31                                                                                              
           echo $thin_b                                                                                               
           echo "$(printf '%*s' $(($center_position-10)) ) Plotting time-average quantities vs time    "     
           echo $thin_b                                                                                               
           tput setaf 0                                                                                               
           python $code_path"plotting_1D_tave_vs_theta.py" $case $start $end $radius $plot_1D_tave_q1 $plot_1D_tave_q2 $plot_1D_tave_q3 $plot_1D_tave_q4 vs_t $plot_1D_tave_save_t $plot_1D_tave_succ
        fi        



    done
done


echo ""                                                                                                            
echo ""                                                                                                                                                                    
echo $(tput setaf 6)$thick_b                                   
echo "$(figlet -f standard -w $((terminal_width-$font_size)) "4.  All Done!!!")"      
echo $(tput setaf 6)$thick_b                                                                                                                                               
sleep 1  

