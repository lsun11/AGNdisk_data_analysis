# load the input file
source input.sh

# Print the banner (require to install figlet)
echo $(tput setaf 1)$thick_b 
echo $(tput setaf 1)$thick_b
#### Somehow the figlet in eddington doesn't have font univers so I use block here...
#echo "$(tput setaf 1)$(tput bold)$(printf '%*s' $center_position)$(figlet -f univers -w $((terminal_width-$font_size)) "${banner_text}")$(printf '%*s' $center_position)$(tput sgr0)"
echo "$(tput setaf 1)$(tput bold)$(printf '%*s' $center_position)$(figlet -f block -w $((terminal_width-$font_size)) "${banner_text}")$(printf '%*s' $center_position)$(tput sgr0)" 
echo "$(tput setaf 1)$(tput bold)$(printf '%*s' $((terminal_width-${#author_text})))${author_text}$(printf '%*s' ${#author_text})$(tput sgr0)"
echo $(tput setaf 1)$thick_b
echo $(tput setaf 1)$thick_b 
echo ""
# set the color to black
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
    for radius in ${radii[@]}; do
	echo $thin_b
	echo $wave_b_short" radius=$radius "$wave_b_short
        echo $thin_b
	echo $start $end $radius
    done
    done




