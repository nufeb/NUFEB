import random
import argparse
import numpy as np
import pickle

parser = argparse.ArgumentParser(description='Create atom definition files')
parser.add_argument('--n', dest='num', action='store',
                   default=1,
                   help='Create atom definition files for NUFEB with --n #files desired (default is 1)')
parser.add_argument('--co2', dest='co2', action='store',
                   default=1e3,
                   help='Set CO2 concentration (mM)')

args = parser.parse_args()
Ks = np.logspace(-2,0,args.num)
CO2MW = 44.01
for n in range(1,int(args.num)+1):
    # atomType = 'cyano'
    n_cells = int(random.uniform(1,100))
    K_co2 = Ks[n-1]#random.uniform(5e0,1e-5)
    min_size = 1.37e-6
    max_size = 1.94e-6
    dimensions = [1e-4,1e-4,1e-5]#x,y,z in meters
    mu = round(0.06/3600,7) #0.047-0.087/hr from Brodderick et al 2019 PCC7942
    InitialConditions = {'CellType': 'cyano', 'StartingCells' : n_cells,'GrowthRate' : mu, 
                         'Nutrients' : {'Concentration' :  {'sub' : 1e-1,'o2' : 9e-3, 'suc' : 1e-20, 'co2' : args.co2*CO2MW*1e-3,'gco2' : 2e-2},
                'State' : {'sub' : 'g','o2' : 'l', 'suc' : 'l', 'co2' : 'l','gco2' : 'g'},
                'xbc' : {'sub' : 'nn','o2' : 'nn', 'suc' : 'nn', 'co2' : 'nn','gco2' : 'pp'},
                'ybc' : {'sub' : 'nn','o2' : 'nn', 'suc' : 'nn', 'co2' : 'nn','gco2' : 'pp'},
                'zbc' : {'sub' : 'nn','o2' : 'nn', 'suc' : 'nn', 'co2' : 'nn','gco2' : 'pp'}},
                'Diff_c' : {'sub' : 0,'o2' : 2.30e-9, 'suc' : 5.2e-10,'co2' : 1.9e-09,'gco2' : 0},
                'K_s' : {'sub' : 3.5e-4,'o2' : 2e-4, 'suc' : 1e-2,'co2' : 1.38e-4,'gco2' : 0},
                'GrowthParams' : {'Yield' : 0.55,'Maintenance' : 0,'Decay' : 0}
                }


    NutesNum = len(InitialConditions['Nutrients']['Concentration'])
    
    
     #yield = 0.55
    
    
    L = [' NUFEB Simulation\r\n\n',f'     {n_cells} atoms \n',
         '     1 atom types \n',f'     {NutesNum} nutrients \n\n',
         f'  0.0e-4   {dimensions[0] :.2e}  xlo xhi \n',f'  0.0e-4   {dimensions[1] :.2e}  ylo yhi \n',
         f'  0.0e-4   {dimensions[2] :.2e}  zlo zhi \n\n', ' Atoms \n\n'
         ]
    
    for i in range(1,n_cells+1):
        size = random.uniform(min_size, max_size)
        x = random.uniform(0+size,dimensions[0]-size)
        y = random.uniform(0+size,dimensions[1]-size)
        z = random.uniform(0+size,dimensions[2]-size)
        L.append(f'     %d 1 {size :.2e}  357 {x :.2e} {y :.2e} {z :.2e} {size :.2e} \n'% (i))
    L.append('\n')
    L.append(' Nutrients \n\n')
    for i,nute in enumerate(InitialConditions['Nutrients']['Concentration'].keys()):
        L.append(f'     %d {nute} {InitialConditions["Nutrients"]["State"][nute]} {InitialConditions["Nutrients"]["xbc"][nute]} {InitialConditions["Nutrients"]["ybc"][nute]} {InitialConditions["Nutrients"]["zbc"][nute]} {InitialConditions["Nutrients"]["Concentration"][nute] :.2e} {InitialConditions["Nutrients"]["Concentration"][nute] :.2e} \n'% (i+1))
       
    L.append('\n')
    L.append(' Type Name \n\n')
    L.append(f'     1 {InitialConditions["CellType"]}  \n\n')
    L.append(' Diffusion Coeffs \n\n')
    for key in InitialConditions['Diff_c'].keys():
        L.append(f'     {key} {InitialConditions["Diff_c"][key]} \n')
    L.append('\n')
    L.append(' Growth Rate \n\n')
    L.append(f'     cyano {InitialConditions["GrowthRate"]} \n\n')
    L.append(' Ks \n\n')    
    k = '     cyano'
    for key in InitialConditions['K_s'].keys():
        k = k + ' ' + str(InitialConditions['K_s'][key])
    k = k + f' \n\n'
    L.append(k) 
    for key in InitialConditions['GrowthParams'].keys():
        L.append(' ' + key + f' \n\n')
        L.append(f'     {InitialConditions["CellType"]} {InitialConditions["GrowthParams"][key]} \n\n')
    f= open(f"atom_{n}.in","w+")
    f.writelines(L)
    dumpfile = open(f"run_{n}.pkl",'wb')
    pickle.dump(InitialConditions,dumpfile)
    dumpfile.close()
