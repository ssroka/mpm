#!/bin/python

import subprocess
import os
import fileinput

# the location of this script
base_dir = os.path.dirname(os.path.realpath(__file__))

com_print_flag = "% print the error between the calculated microphysical quantity and the\n% quantity stated in the corresponding publication"

# will change
RH = 94.0
DT = 3.0
U10 = 20.0
T_s_0 = 27.0
T_a = T_s_0 - DT

# might change
r_0_start = 2000
r_0_inc = 1
r_0_end = 2000

t_final = 2050


# please don't change
p0 = 100000 # Pa
S0 = 34 # ppt
maxit = 1000

com_mp_points_Fig11 = "% micro-physical endpoints from caption in Fig 11 of Andreas 2005"

com_NR = "% the maxium number of iterations before aborting the iterative scheme"

com_print_flag = "% print the error between the calculated microphysical quantity and the\n% quantity stated in the corresponding publication"

# definitely don't change
T_eq_exact = 17.07; # deg C
r_eq_exact = 61.44; # microns
tau_T_exact = 0.0176; # s
tau_r_exact = 303.0; # s

for i in range(200,201):	
	f = open("Exp{}_in.m".format(i),"w")
	f.write("\nsaveDir = 'Exp{}';\n".format(i))

	f.write("\n{}".format(com_print_flag))
	f.write("\nprint_flag = true;\n")

	f.write("\nNayar_flag = true;\n")

	f.write("\n% SSGF String")
	f.write("\nSSGF_str = 'singleDrop';\n")

	f.write("\nRH = {};\n".format(RH))

	f.write("\n% ambient air temperature")
	f.write("\nT_a = {};\n".format(T_a))

	f.write("\n% initial droplet temperature (and local SST)")
	f.write("\nT_s_0 = {};\n".format(T_s_0))

	f.write("\n% initial salinity of drop") 
	f.write("\nS0 = {}; % ppt \n".format(S0))

	f.write("\n% pressure of air") 
	f.write("\np0 = {}; % Pa\n".format(p0))

	f.write("\n% initial droplet radius")
	f.write("\nr_0_vec = [{}:{}:{}]*1e-6; % m\n".format(r_0_start,r_0_inc,r_0_end))

	f.write("\n% 10-m wind speed")
	f.write("\nU10 = {}; % m/s\n".format(U10))

	f.write("\n% integration time")
	f.write("\nt_final = {}; % s\n".format(t_final))

	f.write("\n% Newton-Raphson for r_eq")
	f.write("\nmaxIt = {}; {}\n".format(maxit,com_NR))

	f.write("\nmaxEr_req = max(r_0_vec)*0.00001; % the max error in r_eq")
	f.write("\nmaxEr_uf  = calcVterm_sphere(max(r_0_vec))*0.0001;  % the max error in u_f")
	f.write("\nmaxEr_s   = S0/1000*0.01; % the max error\n") 

	f.write("\n {}".format(com_mp_points_Fig11)) 
	f.write("\nT_eq_exact = {}; % deg C".format(T_eq_exact))
	f.write("\nr_eq_exact = {}* 10^-6; % m".format(r_eq_exact))
	f.write("\ntau_r_exact = {}; % s".format(tau_r_exact))
	f.write("\ntau_T_exact = {}; % s".format(tau_T_exact))


	f.close()









