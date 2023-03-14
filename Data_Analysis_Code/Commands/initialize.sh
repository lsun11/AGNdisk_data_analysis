Directory_plot="/home/lunan/mnt/AGNWedge8_product/PLOTS/"
Directory_checkpoint="/home/lunan/mnt/AGNWedge8_product/Checkpoint/"

declare -a Cases=("Wedge8" "Wedge9B" "Wedge10") 

if [ ! -d "$Directory_plot" ];
then 
	echo -e "Directory for saved plots does not exist! Create one!\n"        
        for val in ${Cases[@]}; do             
            echo "Creating $Directory_plot$val"
            mkdir $Directory_plot$val                                                                                             
        done 
        echo -e "\n\n"
fi

if [ -d "$Directory_plot" ];  
then  
         echo -e "Directory for saved plots exists! Items inside:\n"
         for val in ${Cases[@]}; do
             echo "$Directory_plot$val:"
             ls $Directory_plot$val
         done
         echo -e "\n\n"
fi

if [ ! -d "$Directory_checkpoint" ];   
then                                                                                                                
        echo -e "Directory for saved plots does not exist! Create one!\n"                                           
        for val in ${Cases[@]}; do                                                                             
            echo "Creating $Directory_checkpoint$val"
            mkdir $Directory_checkpoint$val                                                                               
        done
        echo -e "\n\n"
fi                                                                                                                  
                                                                                                                    
if [ -d "$Directory_checkpoint" ];                                                                                        
then                                                                                                                
         echo -e "Directory for saved plots exists! Items inside:\n"                                                
         for val in ${Cases[@]}; do                                                                                 
             echo "$Directory_checkpoint$val:"           
             ls $Directory_checkpoint$val                                                                                      
         done                                                                                                       
         echo -e "\n\n"
fi 