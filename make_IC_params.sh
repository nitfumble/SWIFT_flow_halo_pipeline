#!/bin/bash
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "Making the parameter file for ICICLE.
"

if [[ -e parameterfile.txt ]]; then
	echo "File already exists!
Do you want to make a new file? (y/n)"
	read varname
	if [ $varname != "y" ]; then
		echo "Using the old parameter file. (n = $2)" >&2
	exit 2
	fi
else
	echo "Making a new parameter file."
	rm -f parameterfile.txt
fi

Total_mass=$(( 10 ** $3)) #120 #[10^10 M_Sun]
num_particles=$2
mass=$(awk "BEGIN {printf \"%.5f\",$Total_mass/$num_particles}")

if [ $3 == 2 ]; then
   echo "Model $1
G 4.51846772E-29
m $mass
n $num_particles

#NFW
r_s 26
r_cut 360
truncate True

#NFWX
r_sX 26
r_vir 240.0
d 12

#Hernquist 
a 18

#Einasto
alpha 0.16088252900385014
r2 26

#King
P0 6
r_t 250" > parameterfile.txt 
	echo "The parameterfile is saved"
	echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
elif [ $3 == 3 ]; then
   echo "Model $1
G 4.51846772E-29
m $mass
n $num_particles

#NFW
r_s 70
r_cut 775
truncate True

#NFWX
r_sX 70
r_vir 520.0
d 2

#Hernquist 
a 40

#Einasto
alpha 0.16718100650435025
r2 69.6866

#King
P0 6
r_t 520" > parameterfile.txt 
echo "The parameterfile is saved"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
else
	echo "Model $1
G 4.51846772E-29
m $mass
n $num_particles

#NFW
r_s 220
r_cut 1670
truncate True

#NFWX
r_sX 220
r_vir 1115
d 2

#Hernquist 
a 105

#Einasto
alpha 0.18412102491714508
r2 220

#King
P0 6
r_t 1400.0" > parameterfile.txt 
echo "The parameterfile is saved"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
fi


# UNITS WE WILL USE:
# Time [seconds]; Length [km]; Mass [Solar Masses];
# Alternativly: Time [seconds]; Length [pc]; Mass [kg];
# G: 1.327*10**11 [km^3 s^-2 M⊙^-1] or 4.51667*10**-30 [pc^3 s^-2 M⊙^-1] or 6.67*10**-23 [km^3 s^-2 kg^-1]
# G: 4.5416*10**15 [pc^3 Gyr^-2 M⊙^-1]

# G: 4.51846772E-29 [kpc^3 s^-2 (10^10)M⊙^-1] 
# G: 


# m: 10^10 	[Solar Masses] 
# l: 1 		[pc]
# t: 1 		[Gyr]
# 
# r_s = r200/c; c ≡ Rvir/rs; 
#
# c_MW = ~1
# v ~200km/s
# Rvir 200kpc
# M_MW ~ 10^12 solar masses
