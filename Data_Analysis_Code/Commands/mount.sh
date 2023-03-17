source input.sh

for case in ${cases[@]} ; do
    echo $home_path"DATA"$case
    if [[ $(find $home_path"DATA"$case -mindepth 1 -print -quit) ]]; then
        echo "Directory $home_pathDATA$case is already mounted."
    else
        echo $mid_b
        echo "Connecting Data for case: $case"
        echo $username@128.111.9.179:"/data2/AGN"$case"/History" $home_path"DATA"$case  
        sshfs $username@128.111.9.179:"/data2/AGN"$case"/History" $home_path"DATA"$case 
        echo  "Connection of Data for case: $case successful!!"
        echo $mid_b
    fi
done    
