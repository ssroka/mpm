#!/bin/python

import subprocess
import os
import fileinput

# the location of this script
base_dir = os.path.dirname(os.path.realpath(__file__))

# execute createCase
#runs = [1]
#runs.extend(list(range(2,7)))
#for i in runs:
#    shscript = subprocess.Popen([base_dir+"/createNewCase.sh","bathtub_20_{}0_fine_F".format(i),"bathtub_fine_{}0_U10".format(i)],stdin=subprocess.PIPE)

##  copyfiles
#sn1 = "../simple_coarse_{}0".format(6)
#sn1 = "../bathtub_20_{}0_fine_F".format(1)
#sn1 = "../bathtub_fine_{}0_U10".format(1)
##"/system/controlDict"
# "/system/setFieldsDict"
# "/constant/phaseProperties","/constant/thermophysicalProperties.liquid","/constant/thermophysicalProperties.gas"
# "/top6.setSet"
# "/Allrun"
# "/system/fvOptions"
# "/system/defineGlobalVars"
# "/system/setFieldsDict"
# "/system/blockMeshDict"
#sn1 = "../simple_fine_{}0".format(6)
#files = ["/system/controlDict"]
#for i in range(2,7):
#    #sn1 = "../simple_coarse_{}0".format(i)
#    #sn2 = "simple_coarse_{}0".format(i)
#    #sn2 = "bathtub_20_{}0_coarse_F".format(i)
#    #sn2 = "simple_fine_{}0".format(i)
#    #sn2 = "bathtub_20_{}0_fine_F".format(i)
#    sn2 = "bathtub_fine_{}0_U10".format(i)
#    os.chdir(sn2)
#    for j in range(len(files)):
#        os.system("echo \"Yes\" | cp "+sn1+files[j]+" "+"."+files[j])
#    os.chdir(base_dir)

## create old_log dirs in all Case dirs
#for i in range(1,7):
#    #sn = "bathtub_20_{}0_coarse_F".format(i)
#    #sn = "simple_coarse_{}0".format(i)
#    sn = "simple_fine_{}0".format(i)
#    os.chdir(sn)
#    if not os.path.exists("old_log"):
#        os.makedirs("old_log")
#    os.chdir(base_dir)

## rm files in old_log
#for i in range(1,7):
#    #sn2 = "bathtub_20_{}0_coarse_F".format(i)
#    #sn2 = "simple_fine_{}0".format(i)
#    sn2 = "simple_coarse_{}0".format(i)
#    os.chdir(sn2+"/old_log")
#    os.system("ls -lt")
#    os.system("rm log.icoReactingMultiphaseInterFoam[0-9]")
#    os.system("ls -lt")
#    os.chdir(base_dir)

## rm files in VTK
#for i in range(1,7):
#    #sn2 = "bathtub_20_{}0_coarse_F".format(i)
#    sn2 = "simple_fine_{}0".format(i)
#    os.chdir(sn2+"/VTK")
#    os.system("ls -lt")
#    os.system("rm *.vtk")
#    os.system("ls -lt")
#    os.chdir(base_dir)

## coarse copy log.ico... into old_log dir
#sn1 = "log.icoReactingMultiphaseInterFoam"
#N = "1"
#runs = [1]
#runs.extend(list(range(2,7)))
#for i in runs:
#    #sn2 = "bathtub_20_{}0_coarse_F".format(i)
#    #sn2 = "simple_coarse_{}0".format(i)
#    sn2 = "simple_fine_{}0".format(i)
#    os.chdir(sn2)
#    os.system("mv "+sn1+" old_log/"+sn1+N)
#    os.chdir(base_dir)

##  fine copy log.ico... into old_log dir
#sn1 = "log.icoReactingMultiphaseInterFoam"
#N = "3"
#runs = [1]
#runs.extend(list(range(2,7)))
#for i in runs:
#    sn2 = "bathtub_20_{}0_fine".format(i)
#    os.chdir(sn2)
#    os.system("mv "+sn1+" old_log/"+sn1+N)
#    os.chdir(base_dir)

## coarse execute foamToVTK
for i in range(1,7):
    #sn = "simple_coarse_{}0".format(i)
    #sn = "simple_coarse_{}0_fast".format(i)
    sn = "simple_fine_{}0".format(i)
    #sn = "bathtub_20_{}0_coarse_F".format(i)
    #sn = "bathtub_20_{}0_fine_F".format(i)
    #sn = "bathtub_fine_{}0_U10".format(i)
    os.chdir(sn)
    shscript = subprocess.run(["foamToVTK","-latestTime"],stdin=subprocess.PIPE)
    os.chdir(base_dir)

## fine execute foamToVTK
#for i in range(1,7):
#    #sn = "bathtub_20_{}0_fine".format(i)
#    sn = "simple_fine_{}0".format(i)
#    os.chdir(sn)
#    #shscript = subprocess.run(["foamToVTK","-latestTime"],stdin=subprocess.PIPE)
#    shscript = subprocess.run(["foamToVTK"],stdin=subprocess.PIPE)
#    os.chdir(base_dir)

## execute Allclean coarse
#for i in range(2,7):
#    #sn = "simple_fine_{}0".format(i)
#    #sn = "simple_coarse_{}0".format(i)
#    sn = "bathtub_20_{}0_fine_F".format(i)
#    os.chdir(sn)
#    shscript = subprocess.Popen([base_dir+"/"+sn+"/Allclean"],stdin=subprocess.PIPE)
#    shscript.wait()
#    os.chdir(base_dir)

## execute Allclean fine
#for i in range(1,7):
#    sn = "bathtub_20_{}0_fine".format(i)
#    os.chdir(sn)
#    shscript = subprocess.Popen([base_dir+"/"+sn+"/Allclean"],stdin=subprocess.PIPE)
#    shscript.wait()
#    os.chdir(base_dir)

## simple  launch all the .run scripts with the same "run description" on the same No of cores 
#runs = [1]
#runs.extend(list(range(2,7)))
#for i in runs:
#    #sn = "simple_fine_{}0".format(i)
#    #sn = "simple_coarse_{}0".format(i)
#    #sn = "bathtub_20_{}0_fine_F".format(i)
#    #sn = "simple_coarse_{}0_fast".format(i)
#    sn = "bathtub_fine_{}0_U10".format(i)
#    os.chdir(sn)
#    fn =  "/"+sn+"/"+sn+".run"
#    shscript = subprocess.Popen([base_dir+fn,"try U10 sim","1"],stdin=subprocess.PIPE)
#    shscript.wait()
#    os.chdir(base_dir)

## coarse_F launch all the .run scripts with the same "run description" on the same No of cores 
#runs = [1]
#runs.extend(list(range(2,7)))
#for i in runs:
#    sn = "bathtub_20_{}0_coarse_F".format(i)
#    os.chdir(sn)
#    fn =  "/"+sn+"/"+sn+".run"
#    shscript = subprocess.Popen([base_dir+fn,"try new phase properties","1"],stdin=subprocess.PIPE)
#    shscript.wait()
#    os.chdir(base_dir)

## RESTART coarse launch all the .run scripts with the same "run description" on the same No of cores 
#runs = [1]
#runs.extend(list(range(2,7)))
#for i in runs:
#    #sn = "bathtub_20_{}0_coarse_F".format(i)
#    #sn = "simple_coarse_{}0".format(i)
#    sn = "simple_fine_{}0".format(i)
#    os.chdir(sn)
#    fn =  "/"+sn+"/"+sn+"_RESTART.run"
#    shscript = subprocess.Popen([base_dir+fn,"RESTART for UBUA","1"],stdin=subprocess.PIPE)
#    shscript.wait()
#    os.chdir(base_dir)

## RESTART fine launch all the .run scripts with the same "run description" on the same No of cores 
#runs = [1]
#runs.extend(list(range(2,7)))
#for i in runs:
#    sn = "bathtub_20_{}0_fine".format(i)
#    os.chdir(sn)
#    fn =  "/"+sn+"/"+sn+"_RESTART.run"
#    shscript = subprocess.Popen([base_dir+fn,"RESTART new forcing most failed before t=10s","1"],stdin=subprocess.PIPE)
#    shscript.wait()
#    os.chdir(base_dir)

## fine launch all the .run scripts with the same "run description" on the same No of cores 
#runs = [1]
#runs.extend(list(range(3,7)))
#runs = list(range(1,7))
#for i in runs:
#    sn = "bathtub_20_{}0_fine_F".format(i)
#    os.chdir(sn)
#    fn =  "/"+sn+"/"+sn+".run"
#    shscript = subprocess.Popen([base_dir+fn,"try new domain grid spacing and new fvOptions","1"],stdin=subprocess.PIPE)
#    shscript.wait()
#    os.chdir(base_dir)


# ======= change a single line in a file =======
#for i in range(1,7):
#    sn = "bathtub_20_{}0_coarse".format(i)
#    def_file = sn+"/system/defineGlobalVars"
#    print(def_file)
#    with fileinput.FileInput(def_file,inplace=True,backup='.bak') as file:
#        for line in file:
#            print(line.replace("U10  10.0;", "U10  {}0.0;".format(i)),end='')

# ======= change a single line in a file =======
#str_old ="endTime         100.0;"
#str_old = "writeInterval   500;"
#file_to_change = "/system/controlDict"
#for i in range(1,7):
#    sn = "bathtub_20_{}0_coarse".format(i)
#    def_file = sn+file_to_change
#    print(def_file)
#    #str_new ="endTime         500.0;"
#    str_new = "writeInterval   1000;"
#    with fileinput.FileInput(def_file,inplace=True,backup='.bak') as file:
#        for line in file:
#            print(line.replace(str_old,str_new),end='')

# ======= change a single line in a file =======
#str_old ="endTime         500.0;"
#str_old = "writeInterval   500;"
#file_to_change = "/system/fvOptions"
#for i in range(1,7):
#    sn = "bathtub_20_{}0_coarse_F".format(i)
#    def_file = sn+file_to_change
#    print(def_file)
#    str_old = "uSource[cells_[i]]  += 1.0*(vector({}0.0, 0.0, 0.0) - U[cells_[i]]);".format(i)
#    str_new = "uSource[cells_[i]]  += -1.0*(vector({}0.0, 0.0, 0.0) - U[cells_[i]]);".format(i)
#    #str_new ="endTime         1000.0;"
#    #str_new = "writeInterval   1000;"
#    with fileinput.FileInput(def_file,inplace=True,backup='.bak') as file:
#        for line in file:
#            print(line.replace(str_old,str_new),end='')

# ======= change a single line in a file =======
#str_old ="endTime         50.0;"
#str_old = "writeInterval   500;"
#file_to_change = "/system/controlDict"
#for i in range(1,7):
#    sn = "bathtub_20_{}0_fine".format(i)
#    def_file = sn+file_to_change
#    print(def_file)
#    str_new ="endTime         500.0;"
#    str_new = "writeInterval   1000;"
#    with fileinput.FileInput(def_file,inplace=True,backup='.bak') as file:
#        for line in file:
#            print(line.replace(str_old,str_new),end='')

# ======= change a single line in a file =======
#str_old ="endTime         50.0;"
#str_old = "writeInterval   500;"
#str_old = "#SBATCH --time=24:00:00"
#for i in range(1,7):
#    sn = "bathtub_20_{}0_fine".format(i)
#    file_to_change = "/"+sn+"_RESTART.slurm"
#    def_file = sn+file_to_change
#    print(def_file)
#    str_new ="endTime         500.0;"
#    str_new = "writeInterval   1000;"
#    str_new = "#SBATCH --time=48:00:00"
#    with fileinput.FileInput(def_file,inplace=True,backup='.bak') as file:
#        for line in file:
#            print(line.replace(str_old,str_new),end='')

# get the names of the times timesteps

#for i in range(1,2):
##    sn = "bathtub_20_{}0_fine".format(i)
#	sn = "simple_30_-1"
#	os.chdir(sn)
#	allFiles = os.listdir(".")
#	if not os.path.exists("VTK"):
#		os.makedirs("VTK")
#	f = open("./VTK/filenames.txt","w")
#	for j in allFiles:
#		if j[-1].isdigit():
#			if j=="0":
#				f.write(j+",0\n")
#			else:
#				os.chdir(j+"/uniform")
#				with open("time","r") as ft:
#					line = ft.readline()
#					while line:
#						if len(line)>5 and line[0:5]=="value":
#							f.write(line.split("e")[1].strip().split(";")[0]+",")
#						if len(line)>5 and line[0:5]=="index":
#							f.write(line.split("x")[1].strip().split(";")[0]+"\n")
#						line = ft.readline()
#		os.chdir(base_dir+"/"+sn)
#	f.close()
#	os.chdir(base_dir)













