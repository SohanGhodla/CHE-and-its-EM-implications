#!/usr/bin/env python
import numpy as np
import os

mass = np.arange(15,300,5.0)


z = ['zem5', 'zem4', '6zem4', '3zem4', 'z001', 'z002', 'z003', 'z004', 'z005', 'z006', 'z007', 'z008', 'z009', 'z010']
Z = [0.00001, 0.0001, 0.0006, 0.0003, 0.001, 0.002, 0.003, 0.004,  0.005, 0.006, 0.007, 0.008, 0.009, 0.010]
#z =   ['z040', 'z020', 'z010', 'z008', 'z007', 'z006', 'z005', 'z004', 'z003', 'z002', 'z001', 'zem4', 'zem5']
#Z =   [0.04,   0.02,    0.01,   0.008,  0.007,  0.006,  0.005,  0.004,  0.003,  0.002, 0.001, 0.0001, 0.00001]

rotations = np.arange(0.36,0.98,0.02)
print(rotations)

for rot in rotations:
  rotDirname = "Rot-" + str(round(rot,2))
  for loopz in z:
    for m in mass: 
      dirname = './Makemodels_MS/' +  str(rotDirname) + '/' + str(loopz) + '/' + str(loopz) + '-' + str(m)
      if os.path.isdir(dirname):
        if os.path.isdir(dirname + '/LOGS'):
          #print('found')
          listdirs = os.listdir(dirname)
          for a in range(len(listdirs)):
            if listdirs[a] != 'LOGS':
              os.system('rm -r ' + dirname + '/' + str(listdirs[a]))
              print(dirname + '/' + listdirs[a])
        
            # else:
            #   os.system('rm -r ' + str(dirname))
    
    
