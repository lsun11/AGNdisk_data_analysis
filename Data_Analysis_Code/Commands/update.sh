#### This is the most important command file!! It calls all necessary function possible!!!
case=Wedge8
radius=120

if [[ $case == "Wedge8" ]]; then
   start=2000
   end=3658
elif [[ $case == "Wedge8_2" ]]; then                                                                                                              
   start=2000                 
   end=3658
elif [[ $case == "Wedge9B" ]]; then                                                                                                              
   start=0 
   end=1157
elif [[ $case == "Wedge9B_2" ]]; then
   start=0                  
   end=1213 
else
   echo "Error"
fi


if [[ $radius = '120' ]]
then
    p_gas=pg_r1
elif [[ $radius = '140' ]]  
then         
    p_gas=pg_r2
elif [[ $radius = '160' ]]
then
    p_gas=pg_r3
elif [[ $radius = '180' ]] 
then
    p_gas=pg_r4
else
    echo "Error"
fi

code_path=/home/lunan/mnt/AGNWedge8_product/Data_Analysis_Code

load_time=True
ph_mode=-1
verbose=True

### Make t-average pressure vs theta plot (PB total / pg / Ek2 / Er) at r = 120, 140, 160, 180
##python $code_path/plotting_1D_tave_vs_theta.py $case $start $end $radius PB_total $p_gas Ek2 Er vs_theta $load_time $verbose

### Make quantity-vs-t-vs-theta plots at specific raidus with/without photosphere
python $code_path/plotting_t_vs_th_vs_quant.py $case rho sigma $radius Y $start $end $ph_mode $load_time $verbose