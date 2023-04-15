#!/bin/bash
NOW=$(date +%d%h%Y-%H%M%S)
#Arguments: Model name; #Particles; ;[T/F] Show results  
echo "-------------------------------------------"
echo "Running the SWIFT Pipeline for a halo using the $1 profile with a mass ejection event.
"
#Checking the arguments input.
if [[ $# -ne 3 ]]; then
	echo "Illegal number of arguments;
Please input [Model name] [#particles] [Plot Results (T/F)]" >&2
	exit 2
fi

case $1 in
    NFW|NFWX|Hernquist|King|Einasto) echo "Model is valid."	;;
    *) echo "$1 is an invalid model. Please use: NFW, NFWX ,Hernquist, King or Einasto.">&2
	exit 2 ;;
esac

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Program is set up, Ready to run!

[ ! -d "SwiftRUN_$NOW" ] && mkdir "SwiftRUN_$NOW"
[ ! -d "SwiftRUN_$NOW/Figures" ] && mkdir "SwiftRUN_$NOW/Figures"
[ ! -d "SwiftRUN_$NOW/Results" ] && mkdir "SwiftRUN_$NOW/Results"
[ ! -d "Final_Results" ] && mkdir "Final_Results"
[ ! -d "Final_Results/$1_$NOW" ] && mkdir "Final_Results/$1_$NOW"
cd SwiftRUN_$NOW

for halo_mass in {2..4}
do
solmas=$((10 + $halo_mass))
echo "Running the simulation for a halo of mass: 10^$solmas [M_Sun]"

#Making the parameter file for ICICLE.
#Outputs parameterfile.txt: File with all the info to make the particle IC's given a model.
bash ../make_IC_params.sh $1 $2 $halo_mass

echo "Generating $1 ICs "
#Making the IC's for the given profiles
python3 ../ICICLE.py parameterfile.txt outputfile.csv



#Read the box size:
params=($(<"Params.txt"))
boxsize=${params[0]}
mass=${params[1]}

#Make the IC.hdf5 file from the outputfile of ICICLE.
python3 ../Make_IC_HDF5_File.py outputfile.csv $1 $boxsize $mass $2

#Relaxation run for the IC's
swift --self-gravity --threads=16 Relax.yml 2>&1 | tee output.log

mv IC_$1.hdf5 OLD_IC.hdf5
mv IC_0001.hdf5 IC_$1.hdf5
rm IC_0000.hdf5
rm IC.xmf
rm dependency_graph_0.csv
rm output.log
rm task_level_0.txt
rm timesteps_16.txt
rm unused_parameters.yml
rm used_parameters.yml
rm statistics.txt


for Ejection_time in 4 5.5 7 
do
echo "Running the simulation with the ejection event happening at $Ejection_time Gyr"

#Make .yml File.
python3 ../create_yml.py isolated_galaxy.yml $1 $boxsize $2 $Ejection_time
#2 5 10
for mass_frac in {0..10}
do

#for ps_bool in True False
#do
ps_bool=True

QQ=${ps_bool:0:1}
run_iter=$solmas-$mass_frac-$Ejection_time-$QQ
echo "Run: $run_iter"
[ ! -d "Run_$run_iter" ] && mkdir "Run_$run_iter"

cp -f isolated_galaxy.yml Run_$run_iter/isolated_galaxy.yml
cp -f isolated_galaxy_ejection.yml Run_$run_iter/isolated_galaxy_ejection.yml
cp -f IC_$1.hdf5 Run_$run_iter/IC_$1.hdf5
cp -f parameterfile.txt Run_$run_iter/parameterfile.txt
cd Run_$run_iter

#Run SWIFT relaxation run:
swift --self-gravity --threads=16 isolated_galaxy.yml 2>&1 | tee output.log

echo "The Ejection Event."

#the Ejection event:
python3 ../../Ejection_event.py $mass_frac $boxsize $ps_bool

#Run SWIFT run after mass "ejection" from the centre:
swift --self-gravity --threads=16 isolated_galaxy_ejection.yml 2>&1 | tee output_eject.log

#Plot the results if plotting is True.
echo "Saving data of the conservation of total angular momentum" 
python3 ../../Plot_Func/angularmomentum.py "../Figures/" $1-$run_iter
#echo "Saving data of the evolution of the profile over time." 
#python3 ../Plot_Func/profiles_gif.py "Figures/" $1
echo "Saving data of the evolution of the positions over time." 
python3 ../../Plot_Func/halo_vision.py $1-$run_iter $3
	
cd ..
#done
done
done
#cd ..

rm -f parameterfile.txt
rm -f Params.txt
rm -f isolated_galaxy.yml
rm -f IC_NFW.hdf5
rm -f outputfile.csv
rm -f isolated_galaxy_ejection.yml

done

python3 ../Plot_Func/final_plot.py $1_$NOW  

echo "End of the SWIFT Pipeline"
echo "-------------------------------------------"
