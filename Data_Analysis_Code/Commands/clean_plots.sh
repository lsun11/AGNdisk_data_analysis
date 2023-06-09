Data_dir=$(dirname $(dirname $(pwd)))"/"
Plot="PLOTS/"

Plot_dir=$Data_dir$Plot                                                                   
Checkpoint_dir=$Data_dir$Checkpoint 
echo $Plot_dir

read -p "Delete all plots in "$Plot_dir", are you sure? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    rm -rf $Plot_dir/*/*
fi
echo "All plots in "$Plot_dir" deleted."