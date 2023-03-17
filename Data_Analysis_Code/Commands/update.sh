# load the input file
source input.sh

# Print the banner (require to install figlet)
echo $(tput setaf 1)$thick_b 
echo $(tput setaf 1)$thick_b
echo $(tput setaf 1)$thick_b 
#### Somehow the figlet in eddington doesn't have font univers so I use block here...
#echo "$(tput setaf 1)$(tput bold)$(printf '%*s' $center_position)$(figlet -f univers -w $((terminal_width-$font_size)) "${banner_text}")$(printf '%*s' $center_position)$(tput sgr0)"
echo "$(tput setaf 1)$(tput bold)$(printf '%*s' $center_position)$(figlet -f block -w $((terminal_width-$font_size)) "${banner_text}")$(printf '%*s' $center_position)$(tput sgr0)" 
echo "$(tput setaf 1)$(tput bold)$(printf '%*s' $((terminal_width-${#author_text}-17)))${author_text}$(printf '%*s' ${#author_text} )$(tput sgr0)"
echo "$(tput setaf 1)$(tput bold)$(printf '%*s' $((terminal_width-${#current_time}-15)))"Running Time: "${current_time}$(printf '%*s' ${#current_time} )$(tput sgr0)"
echo $(tput setaf 1)$thick_b
echo $(tput setaf 1)$thick_b
echo $(tput setaf 1)$thick_b  
echo ""
echo "" 
sleep 2

echo $(tput setaf 6)$thick_b                                                                                                                                                                                                  
echo "$(figlet -f standard -w $((terminal_width-$font_size)) "Initialize!!!")"                                                                                                                                        
echo $(tput setaf 6)$thick_b                                                                    
tput setaf 0 
./initialize.sh
sleep 1

echo ""
echo ""
echo $(tput setaf 6)$thick_b                                                                                                                                                                                                                   
echo "$(figlet -f standard -w $((terminal_width-$font_size)) "Mount Data!!!")" 
echo $(tput setaf 6)$thick_b 
tput setaf 0
./mount.sh

echo ""
echo ""
echo $(tput setaf 6)$thick_b 
echo "$(figlet -f standard -w $((terminal_width-$font_size)) "Start Making Plots!!!")"
echo $(tput setaf 6)$thick_b                                                                                                                                                                                                  
sleep 1

tput setaf 0
for case in ${cases[@]} ; do
    echo ""
    if [[ $case == "Wedge8" ]]; then  
        echo $mid_b
	echo $vert_b_short" Case Wedge8 "$vert_b_short
        echo $mid_b
	start=2000
	end=3658
    elif [[ $case == "Wedge8_2" ]]; then
      	echo $mid_b
        echo $vert_b_short" Case Wedge8_2 "$vert_b_short
	echo $mid_b
	start=2000
	end=3658
    elif [[ $case == "Wedge9B" ]]; then
        echo $mid_b
        echo $vert_b_short" Case Wedge9B "$vert_b_short
        echo $mid_b
	start=0
	end=1157
    elif [[ $case == "Wedge9B_2" ]]; then
        echo $mid_b
        echo $vert_b_short" Case Wedge9B_2 "$vert_b_short
        echo $mid_b
	start=0
	end=1213
    else
	echo "Error"
	exit 1;
    fi

    if [[ $plot_sin_cos_quant =~ ^[Yy]$ ]]; 
    then
        echo $mid_b
        echo "$(printf '%*s' $(($center_position-22)) ) Plotting poloidal distribution plot of azimuthally-averaged quantity! "
        echo $mid_b
        for iter in ${plot_sin_cos_iter[@]}; do
            python $code_path"plotting_costh_vs_sinth_vs_quant.py" $case $plot_sin_cos_q1 $plot_sin_cos_q2 $iter $plot_sin_cos_plot_phot $plot_sin_cos_phot_mode $plot_sin_cos_succ;
        done
    fi
    

    for radius in ${radii[@]}; do	
        echo $thin_b
	echo $wave_b_short" radius=$radius "$wave_b_short
        echo $thin_b                 
        if [[ $plot_t_th_dens =~ ^[Yy]$ ]];
        then
           echo "$(printf '%*s' $(($center_position-22)) ) Plotting density vs times vs theta (with photosphere)   "         
           echo $mid_b 
           python $code_path"plotting_t_vs_th_vs_quant.py" $case rho $plot_t_th_sigma None $radius Y $start $end $plot_t_th_phot_mode $plot_t_th_phot_save_t $plot_t_th_phot_succ;
        fi
    done
done


echo ""                                                                                                            
echo ""                                                                                                                                                                    
echo $(tput setaf 6)$thick_b                                   
echo "$(figlet -f standard -w $((terminal_width-$font_size)) "All Done!!!")"      
echo $(tput setaf 6)$thick_b                                                                                                                                               
sleep 1  

