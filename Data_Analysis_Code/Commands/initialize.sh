source input.sh

echo $thick_b                                                                                                                                                                                                       
echo "  Let's check if plots and checkpoints directories are in correct positions!"                                                                                                                                                
echo $thick_b 
echo ""

### Check plot path ###
if [ -d "$plot_path" ];                                                                                                                                                                                             
then                                                                                                                                                                                                                
         echo $(tput setaf 6)$wave_b                                                                                                                                                                               
         echo -e "  Directory for saved plots exists! Items inside:"                                                                                                                                                
         echo $(tput setaf 6)$wave_b                                                                                                                                                                               
         for val in ${cases[@]}; do                                                                                                                                                                                 
             echo $(tput setaf 4)$thin_b                                                                                                                                                                            
             echo "$plot_path$val:"     
             echo $(tput setaf 4)$thin_b                                                                                                                                                                            
             ls $plot_path$val                                                                                                                                                                                      
         done                                                                                                                                                                                                       
         echo -e "\n" 
fi 

if [ ! -d "$plot_path" ];
then 
        echo $(tput setaf 6)$wave_b 
	echo -e "  Directory for saved plots does NOT exist! Create one!"        
        echo $(tput setaf 6)$wave_b                                                                                                                                                                                        
        mkdir $plot_path
        for val in ${cases[@]}; do             
            echo $(tput setaf 4)$thin_b
            echo "  Creating $plot_path$val"
            echo $(tput setaf 4)$thin_b 
            mkdir $plot_path$val                                                                                             
        done 
        echo -e "\n"
fi


### Check checkpoint path ###
if [ -d "$checkpoint_path" ];                                                                                                                                                                                       
then                                                                                                                                                                                                                
         echo $(tput setaf 6)$wave_b                                                                                                                                                                               
         echo -e "  Directory for checkpoints exists! Items inside:"                                                                                                                                                
         echo $(tput setaf 6)$wave_b                                                                                                                                              
         for val in ${cases[@]}; do                                                                                                                                                                                 
             echo $(tput setaf 4)$thin_b                                                                                                                                                                            
             echo "$checkpoint_path$val:"                                                                                                                                                                           
             echo $(tput setaf 4)$thin_b                                                                                                                                                                            
             ls $checkpoint_path$val                                                                                                                                                                                
         done                                                                                                                                                                                                       
         echo -e "\n"
fi

if [ ! -d "$checkpoint_path" ];   
then        
        echo $(tput setaf 6)$wave_b                                                                                                         
        echo -e "  Directory for checkpoints does NOT exist! Create one!"                                           
        echo $(tput setaf 6)$wave_b 
        mkdir $checkpoint_path
        for val in ${cases[@]}; do                                                                             
            echo $(tput setaf 4)$thin_b 
            echo "  Creating $checkpoint_path$val"
            echo $(tput setaf 4)$thin_b 
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

                                                                                         
tput setaf 0
echo $thick_b
echo "  Plots and checkpoints directories are now in correct positions!" 
echo $thick_b 
