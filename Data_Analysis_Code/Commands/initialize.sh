source input.sh

tput setaf 6
echo $mid_b                                                                                                                       
echo "  Let's check if plots and checkpoints directories are in correct positions!"                                                                                                                                                
echo $mid_b 


### Check plot path ###
if [ -d "$plot_path" ];                                                                                                                                                                                             
then                                                                                                                                                                                                                
         echo $(tput setaf 24)$thin_b         
         echo -e "  Directory for saved plots exists! Items inside:"                                                                                                                                                
         echo $(tput setaf 24)$thin_b                                                                                                                                                                               
         for val in ${cases[@]}; do                                                                                                                                                                                 
             echo $(tput setaf 39)$wave_b_short
             echo "$plot_path$val:"     
             echo $(tput setaf 39)$wave_b_short
             tput setaf 0
             ls $plot_path$val                                                                                                                                                                                      
         done                                                                                                                                                                                                       
         echo -e "\n" 
fi 

if [ ! -d "$plot_path" ];
then 
        echo $(tput setaf 24)$thin_b 
	echo -e "  Directory for saved plots does NOT exist! Create one!"        
        echo $(tput setaf 24)$thin_b                                                                                                                                                                                        
        mkdir $plot_path
        for val in ${cases[@]}; do             
            echo $(tput setaf 39)$wave_b_short
            echo "  Creating $plot_path$val"
            echo $(tput setaf 39)$wave_b_short 
            mkdir $plot_path$val                                                                                             
        done 
        echo -e "\n"
fi


### Check checkpoint path ###
if [ -d "$checkpoint_path" ];                                                                                                                                                                                       
then                                                                                                                                                                                                                
         echo $(tput setaf 24)$thin_b 
         echo -e "  Directory for checkpoints exists! Items inside:"                                                                                                                                                
         echo $(tput setaf 24)$thin_b                                                                                                                                              
         for val in ${cases[@]}; do                                                                                                                                                                                 
             echo $(tput setaf 39)$wave_b_short
             echo "$checkpoint_path$val:"                                                                                                                                                                           
             echo $(tput setaf 39)$wave_b_short                                                                                                                                                                      
             tput setaf 0
             ls $checkpoint_path$val                                                                                                                                                                                
         done                   
         echo -e "\n"
fi

if [ ! -d "$checkpoint_path" ];   
then        
        echo $(tput setaf 24)$thin_b                                                                                                         
        echo -e "  Directory for checkpoints does NOT exist! Create one!"                                           
        echo $(tput setaf 24)$thin_b 
        mkdir $checkpoint_path
        for val in ${cases[@]}; do                                                                             
            echo $(tput setaf 39)$wave_b 
            echo "  Creating $checkpoint_path$val"
            echo $(tput setaf 39)$wave_b 
            mkdir $checkpoint_path$val                                                                               
        done
        echo -e "\n"
fi          



### Check Data path ###
for val in ${cases[@]}; do 
    echo $val
    if [ ! -d $home_path"DATA"$val ];
    then 
        echo $(tput setaf 6)$wave_b                                                                                                                                                      
        echo -e "  Directory for Case:$val data does NOT exist! Create one!"                                                                                                                
        echo $(tput setaf 6)$wave_b
        mkdir $home_path"DATA"$val 
    fi
done
echo -e "\n"    

                                                                                         
tput setaf 6
echo $mid_b
echo "  Plots and checkpoints directories are now in correct positions!" 
echo $mid_b 
tput setaf 0