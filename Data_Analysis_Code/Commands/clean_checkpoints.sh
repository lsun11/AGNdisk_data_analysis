Data_dir=$(dirname $(dirname $(pwd)))"/" 
Checkpoint="Checkpoint/"
Checkpoint_dir=$Data_dir$Checkpoint 
echo $Checkpoint_dir

read -p "Delete all plots in "$Checkpoint_dir", are you sure? " -n 1 -r  
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    rm -rf $Checkpoint_dir/*/*
fi
echo "All checkpoint files in "$Checkpoint_dir" deleted."