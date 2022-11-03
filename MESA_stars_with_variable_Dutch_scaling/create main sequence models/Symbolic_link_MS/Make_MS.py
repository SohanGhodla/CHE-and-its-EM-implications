#!/usr/bin/env python
import os
import numpy as np

mass1 = np.arange(15,145,10.0)
mass2 = np.arange(150,226,25.0)
mass = mass1.tolist() + mass2.tolist()

z = ['zem5', 'zem4',  '3zem4', '6zem4', 'z001', 'z002', 'z003', 'z004', 'z005', 'z006']
Z = [0.00001, 0.0001,  0.0003, 0.0006, 0.001, 0.002, 0.003, 0.004, 0.005,  0.006]

rotations = np.arange(0.36,0.98,0.02)

#os.system('rm -r ../Makemodels_MS') #**********************remove Makemodels_MS********************
# Datafile = open('../Run_data.txt' , 'w')
count_symbolic = -1

for rot in rotations:
  count = -1
  rotDirname = "Rot-" + str(np.round(rot,2))
  for loopz in z:
    count += 1
    os.system('mkdir -p ../Makemodels_MS/' + rotDirname + '/' + str(loopz))
    
    for m in mass:
      count_symbolic += 1
      metal = Z[count]
      Dutch_factor = np.sqrt(metal/0.02)
      
      dirname = str(loopz) + '-' + str(m)
     
      print('Z = ' + str(loopz) + ', mass = '  + str(m) +  ', Omega = ' + str(rot) + ', Symlink = ' + str(count_symbolic) + '\n')
      # Datafile.write('Z = ' + str(loopz) + ', mass = '  + str(m) +  ', Omega = ' + str(rot) + ', Symlink = ' +  str(count_symbolic) + '\n')

      if os.path.isdir('../Makemodels_MS/' + str(rotDirname) + '/' + dirname):
        continue
      else:
        os.system('cp -r ../XX-clean_script/ ../Makemodels_MS/' + str(rotDirname) + '/' + str(loopz) + '/' + str(dirname))
        INLIST = open('../XX-clean_script/inlist_var_input', 'r')
        TFILE = open('../Makemodels_MS/' + str(rotDirname) + '/' + str(loopz) + '/' + str(dirname) + '/inlist_var_input', 'w')
        datalines = INLIST.read()
        for line in datalines.split('\n'):
          check = 0
          if 'saved_model_name' in line:
            TFILE.write('\tsaved_model_name = \'zams_' + str(m) + 'Msun.mod\'' + '\n')
            check = 1
          if 'new_omega_div_omega_crit' in line:
            TFILE.write('\tnew_omega_div_omega_crit = ' + str(np.round(rot,2)) + '\n')
            check = 1
          if 'initial_mass =' in line:
            TFILE.write('\tinitial_mass = ' + str(m) + '\n')
            check = 1
          if 'initial_z =' in line:
            TFILE.write('\tinitial_z = ' + str(metal) + '\n')
            check = 1
          if 'Zbase =' in line:
            TFILE.write('\tZbase = ' + str(metal) + '\n')
            check = 1
          if 'Dutch_scaling_factor =' in line:
            TFILE.write('\tDutch_scaling_factor = ' + str(np.round(Dutch_factor,2)) + '\n')
            check = 1
          #************************************************************
          # if 'overshoot_f(1)' in line and m < 60:
          #   TFILE.write('\tovershoot_f(1) = ' + str(L_ov[0]) + '\n')
          #   check = 1
          # if 'overshoot_f(1)' in line and m >= 60 and m < 140:
          #   TFILE.write('\tovershoot_f(1) = ' + str(L_ov[1]) + '\n')
          #   check = 1
          # if 'overshoot_f(1)' in line and m >= 140: 
          #   TFILE.write('\tovershoot_f(1) = ' + str(L_ov[2]) + '\n')
          #   check = 1
           #************************************************************

           #***********************Relaxing the metallicity to desired Z*******************************
          # relax_initial_Z = .true.
          # new_Z = 0.001d0
          # relax_initial_Y = .true.
          # new_Y = 0.242d0

          # if 'new_Z =' in line:
          #   TFILE.write('\tnew_Z = ' + str(metal) + '\n')
          #   check = 1
          # if 'new_Y =' in line:
          #   TFILE.write('\tnew_Y = ' + str(0.24 + 2*metal) + '\n')
            # check = 1
          #********************************************************************************************
          
          if check < 0.5:
            TFILE.write(line + '\n')
        INLIST.close()
        TFILE.close()
    
        os.system('cp -r ./Save_ZAMS_models/' + str(loopz)  + '/zams_' + str(m) + 'Msun.mod ../Makemodels_MS/' + str(rotDirname) + '/' + str(loopz) + '/' + str(dirname))
        
      os.system('ln -s ../Makemodels_MS/' + str(rotDirname) + '/' + str(loopz) + '/' + str(dirname) + ' ' + str(count_symbolic))

