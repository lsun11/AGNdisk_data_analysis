#read -p "Enter a match word: " USER_INPUT
USER_INPUT=hello
if [[ $USER_INPUT == "hello" ]]; then
     echo "Hello to you as well"
fi
case=Wedge9B_2
#radius=120                 
echo $case, [[ $case == "Wedge9B" ]]                                                                              
#echo $radius                                                                                                      

if [[ $case == "Wedge8" ]]; then
   start=2000
   end=3658
elif [[ $case == "Wedge8_2" ]]; then 
   start=2000
   end=3658
elif [[ $case == "Wedge9B" ]]; then
   echo "No"
   start=0
   end=1157
elif [[ $case == "Wedge9B_2" ]]; then
   echo "Yes"
   start=0
   end=1213
fi
echo $start $end


if [[$case=="Wedge8"]]; then  
   start = 2000                                                                                                      
   end = 3658                                                                                                     
elif [[$case=="Wedge8_2"]]; then             
   start = 2000                                                                                                      
   end = 3658                                                                                                     
elif [[$case=="Wedge9B"]];  then                 
   echo "Yes2" 
   start = 0                                                                                                         
   end = 1157                                                                                                     
elif [[$case=="Wedge9B_2"]]; then                  
   start = 0                                                                                                         
   end = 1213                  
else
   echo "error"
fi 
echo $start