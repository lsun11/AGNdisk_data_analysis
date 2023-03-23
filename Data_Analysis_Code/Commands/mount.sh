source input.sh

for case in ${cases[@]} ; do
    echo $(tput setaf $bl)$thin_b
    echo "Check if data is mounted"
    echo $thin_b
    if [[ $(find $home_path"DATA"$case -mindepth 1 -print -quit) ]]; then
        echo $wave_b
        echo "Directory "$home_path"DATA"$case" is ALREADY mounted."
        echo $wave_b
    else
        echo $mid_b
        echo "Data not mounted! Connecting Data for case: $case"
        echo $thin_b 
        echo "Mount: " $username@128.111.9.179:"/data2/AGN"$case"/History" 
        echo "To:     local:"$home_path"DATA"$case 
        echo ""
        echo "Please type the password!!!!" 
        echo $thin_b 
        sshfs $username@128.111.9.179:"/data2/AGN"$case"/History" $home_path"DATA"$case 
        echo  "Connection of Data for case: $case successful!!"
        echo $mid_b
    fi
    tput setaf 0
done    
