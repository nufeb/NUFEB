import random
import argparse
import numpy as np
import pickle
from string import Template

parser = argparse.ArgumentParser(description='Create atom definition files')
parser.add_argument('--n', dest='num', action='store',
                   default=1,
                   help='Create atom definition files for NUFEB with --n #files desired (default is 1)')
parser.add_argument('--r', dest='reps', action='store',
                   default=1,
                   help='Number of replicates')
parser.add_argument('--c',dest='culture_type',action='store',default='co',
                    help='Set culture conditions with --c (co-culture), --ax-c (cyano), --ax-e (e.coli)')
parser.add_argument('--co2', dest='co2', action='store',
                   default=1e3,
                   help='Set CO2 concentration (mM)')
parser.add_argument('--d', dest='dims', action='store',
                   default=[1e-4,1e-4,1e-5],
                   help='Set simulation box dimensions (m)')
parser.add_argument('--t', dest='timesteps', action='store',
                   default=30000,
                   help='Number of timesteps to run')

args = parser.parse_args()

mu_cyanos = round(0.06/3600,7)
mu_ecw = 2.7e-04
CO2MW = 44.01
for n in range(1,int(args.num)+1):
    dimensions = args.dims#x,y,z in meters
    Replicates = args.reps
    SucRatio = round(random.random(),3)
    if args.culture_type == 'co':
        cell_types = ['cyano','ecw']
        n_cyanos = int(random.uniform(1,100))
        n_ecw = int(random.uniform(1,100))
        n_cells = n_cyanos + n_ecw
    elif args.culture_type == 'ax-c':
        cell_types = ['cyano']
        n_cyanos = int(random.uniform(1,100))
        n_ecw = 0
        n_cells = n_cyanos
    elif args.culture_type == 'ax-e':
        cell_types = ['ecw']
        n_ecw = int(random.uniform(1,100))
        n_cyanos=0
        n_cells = n_ecw
    InitialConditions = {'cyano': {'StartingCells' : n_cyanos,'GrowthRate' : mu_cyanos,
           'min_size' : 1.37e-6, 'max_size' : 1.94e-6, 'Density' : 370,
             'K_s' : {'sub' : 3.5e-4,'o2' : 2e-4, 'suc' : 1e-2,'co2' : 1.38e-4,'gco2' : 0},
            'GrowthParams' : {'Yield' : 0.55,'Maintenance' : 0,'Decay' : 0}},
             'ecw': {'StartingCells' : n_ecw,'GrowthRate' : mu_ecw,
           'min_size' : 8.8e-7, 'max_size' : 1.04e-6, 'Density' : 236,
             'K_s' : {'sub' : 0,'o2' : 1e-3, 'suc' : 2.29,'co2' : 5e-2,'gco2' : 0},
            'GrowthParams' : {'Yield' : 0.43,'Maintenance' : 0,'Decay' : 0}},
            'Nutrients' : {'Concentration' :  {'sub' : 1e-1,'o2' : 9e-3, 'suc' : 1e-20, 'co2' : args.co2*CO2MW*1e-3,'gco2' : 2e-2},
            'State' : {'sub' : 'g','o2' : 'l', 'suc' : 'l', 'co2' : 'l','gco2' : 'g'},
            'xbc' : {'sub' : 'nn','o2' : 'nn', 'suc' : 'nn', 'co2' : 'nn','gco2' : 'pp'},
            'ybc' : {'sub' : 'nn','o2' : 'nn', 'suc' : 'nn', 'co2' : 'nn','gco2' : 'pp'},
            'zbc' : {'sub' : 'nn','o2' : 'nn', 'suc' : 'nn', 'co2' : 'nn','gco2' : 'pp'}},
            'Diff_c' : {'sub' : 0,'o2' : 2.30e-9, 'suc' : 5.2e-10,'co2' : 1.9e-09,'gco2' : 0},
            'Dimensions' : dimensions,'SucRatio' : SucRatio,'Replicates' : Replicates

            }
    
    NutesNum = len(InitialConditions['Nutrients']['Concentration'])
    
    L = [' NUFEB Simulation\r\n\n',f'     {n_cells} atoms \n',
         f'     {len(cell_types)} atom types \n',f'     {NutesNum} nutrients \n\n',
         f'  0.0e-4   {dimensions[0] :.2e}  xlo xhi \n',f'  0.0e-4   {dimensions[1] :.2e}  ylo yhi \n',
         f'  0.0e-4   {dimensions[2] :.2e}  zlo zhi \n\n', ' Atoms \n\n'
         ]
    j = 1
    for c, CellType in enumerate(cell_types,start=1):
        for i in range(j,InitialConditions[CellType]['StartingCells']+j):
            size = random.uniform(InitialConditions[CellType]['min_size'], 
                                  InitialConditions[CellType]['max_size'])
            x = random.uniform(0+size,InitialConditions['Dimensions'][0]-size)
            y = random.uniform(0+size,InitialConditions['Dimensions'][1]-size)
            z = random.uniform(0+size,InitialConditions['Dimensions'][2]-size)
            L.append(f'     %d {c} {size :.2e}  {InitialConditions[CellType]["Density"]} {x :.2e} {y :.2e} {z :.2e} {size :.2e} \n'% (i))
            j += 1

    L.append('\n')
    L.append(' Nutrients \n\n')
    for i,nute in enumerate(InitialConditions['Nutrients']['Concentration'].keys()):
        L.append(f'     %d {nute} {InitialConditions["Nutrients"]["State"][nute]} {InitialConditions["Nutrients"]["xbc"][nute]} {InitialConditions["Nutrients"]["ybc"][nute]} {InitialConditions["Nutrients"]["zbc"][nute]} {InitialConditions["Nutrients"]["Concentration"][nute] :.2e} {InitialConditions["Nutrients"]["Concentration"][nute] :.2e} \n'% (i+1))
    
    L.append('\n')
    L.append(' Type Name \n\n')
    for c, CellType in enumerate(cell_types,start=1):
        L.append(f'     {c} {CellType}  \n')
    L.append('\n')
    L.append(' Diffusion Coeffs \n\n')
    for key in InitialConditions['Diff_c'].keys():
        L.append(f'     {key} {InitialConditions["Diff_c"][key]} \n')
    L.append('\n')
    L.append(' Growth Rate \n\n')
    for CellType in cell_types:
        L.append(f'     {CellType} {InitialConditions[CellType]["GrowthRate"]} \n')
    L.append('\n')
    L.append(' Ks \n\n')   
    for CellType in cell_types:
        k = f'     {CellType}'
        for key in InitialConditions[CellType]['K_s'].keys():
            k = k + ' ' + str(InitialConditions[CellType]['K_s'][key])
        k = k + f' \n'
        L.append(k) 
    L.append('\n')
    for key in InitialConditions["cyano"]['GrowthParams'].keys():
        L.append(' ' + key + f' \n\n')
        for CellType in cell_types:
            L.append(f'     {CellType} {InitialConditions[CellType]["GrowthParams"][key]} \n')
        L.append('\n')
        
            
    L.append('\n\n')
    #write atom definition file
    f= open(f"atom_{n}.in","w+")
    f.writelines(L)
    #write initial conditions pickle file
    dumpfile = open(f"run_{n}.pkl",'wb')
    pickle.dump(InitialConditions,dumpfile)
    dumpfile.close()
    #write Inputscript
    #open the file
    filein = open( 'Inputscript-template.txt' )
    #read it
    src = Template( filein.read() )
    #do the substitution
    result = src.safe_substitute({'n' : n, 'SucRatio' : SucRatio,
                                  'Replicates' : Replicates,'Timesteps' : args.timesteps})
    f= open(f"Inputscript_{n}.lammps","w+")
    f.writelines(result)
    #write slurm script
    #open the file
    filein = open( 'slurm-template.txt' )
    #read it
    src = Template( filein.read() )
    #do the substitution
    result = src.safe_substitute({'n' : n, 'job' : f"NUFEB_cyano{n}",'USER' : 'sakkosjo','reps'  : Replicates})
    f= open(f"Inputscript_{n}.slurm","w+")
    f.writelines(result)
