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

echo "Model $1
G 4.51846772E-29
m 0.00014242566
n $2

#NFW
r_s 15
r_cut 200
truncate True

#NFWX
r_sX 15
r_vir 250
d 20

#Hernquist
a 0.8

#King
P0 1.0
r_t 230

#Einasto
r2 1.0
alpha 0.15" > parameterfile.txt 
echo "The parameterfile is saved"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

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
