#!/bin/bash
#Arguments: Model name; #Particles; ;[T/F] Show results  
echo "-------------------------------------------"
echo "Running the SWIFT Pipeline
"
#Checking the arguments input.
if [[ $# -ne 3 ]]; then
	echo "Illegal number of arguments;
Please input [Model name] [#particles] [Plot Results (T/F)]" >&2
	exit 2
fi

case $1 in
    NFW|NFWX|Hernquist|King|Einasto) echo "Model is valid."	;;
    *) echo "Invalid model. Please use: NFW, NFWX ,Hernquist, King, Einasto.">&2
	exit 2 ;;
esac

#Making the parameter file for ICICLE.
bash make_IC_params.sh $1 $2

[ ! -d "SwiftRUN_$1" ] && mkdir "SwiftRUN_$1"
[ ! -d "SwiftRUN_$1/Figures" ] && mkdir "SwiftRUN_$1/Figures"

#Outputs parameterfile.txt: File with all the info to make the particle IC's given a model.
python ICICLE.py parameterfile.txt outputfile.csv

#Read the box size:
params=($(<"Params.txt"))
boxsize=${params[0]}
mass=${params[1]}

#Make the IC.hdf5 file from the outputfile of ICICLE.
python3 Make_IC_HDF5_File.py outputfile.csv $1 $boxsize $mass

##Optional make/change .yml File.
python3 create_yml.py isolated_galaxy.yml $1 $boxsize

#Move files to folder destination:
mv -f Params.txt SwiftRUN_$1/Params.txt
mv -f isolated_galaxy.yml SwiftRUN_$1/isolated_galaxy.yml
mv -f outputfile.csv SwiftRUN_$1/outputfile.csv
mv -f parameterfile.txt SwiftRUN_$1/parameterfile.txt
mv -f IC_$1.hdf5 SwiftRUN_$1/IC_$1.hdf5

cd SwiftRUN_$1
#Run SWIFT:
swift --external-gravity --self-gravity --threads=16 isolated_galaxy.yml 2>&1 | tee output.log

#Plot the results if plotting is True.
if [[ $3 = T ]]; then
	echo "Do you want to plot the conservation of total angular momentum? (y/n)"
	read varname1
	if [ $varname1 == "y" ]; then
		echo "Making plots of the conservation of total angular momentum" 
		python3 ../Plot_Func/angularmomentum.py "Figures/" $1
	fi

	echo "Do you want to plot the evolution of the mass-profile over time (as .GIF)? (y/n)"
	read varname2
	if [ $varname2 == "y" ]; then
		echo "Making plots of the evolution of the profile over time." 
		python3 ../Plot_Func/profiles_gif.py "Figures/" $1
	fi

	echo "Do you want to plot the evolution of the particle positions over time (as .GIF)? (y/n)"
	read varname3
	if [ $varname3 == "y" ]; then
		echo "Making plots of the evolution of the positions over time." 
		python3 ../Plot_Func/halo_vision.py "Figures/" $1
	fi

#	echo "Do you want to plot the change of vertical and radial profile? (y/n)"
#	read varname
#	if [ $varname != "y" ]; then
#		echo "Making plots of the change of vertical and radial profile" 
#		python3 ../Plot_Func/profilefit.py "Figures/" $1
#	fi

#	echo "Do you want to plot the ? (y/n)"
#	read varname
#	if [ $varname != "y" ]; then
#		echo "Making plots of " 
#		python3 ../
#	fi	
	
else
	echo "No results are plotted, see the output .hdf5 files for the results."
fi
echo "End of the SWIFT Pipeline"
echo "-------------------------------------------"
