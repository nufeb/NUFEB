import random
import argparse

parser = argparse.ArgumentParser(description='Create atom definition files')
parser.add_argument('--n', dest='num', action='store',
                   default=1,
                   help='Create atom definition files for NUFEB with --n #files desired (default is 1)')

args = parser.parse_args()

for n in range(1,int(args.num)+1):
    atomType = 'cyano'
    n_cells = int(random.uniform(1,100))
    K_co2 = random.uniform(1e-1,1e-3)
    min_size = 1.37e-6
    max_size = 1.94e-6
    dimensions = [1e-4,1e-4,1e-5]#x,y,z in meters
    growthRate = round(0.047/3600,7) #0.047-0.087/hr from Brodderick et al 2019 PCC7942
    Nutrients = {'sub' : 1e-1,'o2' : 9e-3, 'suc' : 1e-20, 'co2' : 1.2e-1}
    NutesNum = len(Nutrients)
    Diff_c = {'sub' : 0,'o2' : 2.30e-9, 'suc' : 5.2e-10,'co2' : 1.9e-09}
    K_s = {'sub' : 3.5e-4,'o2' : 2e-4, 'suc' : 3.4,'co2' : K_co2}
    Params = {'Yield' : .3,'Maintenance' : 0,'Decay' : 0}
    
    
    L = [' NUFEB Simulation\r\n\n',f'     {n_cells} atoms \n',
         '     1 atom types \n',f'     {NutesNum} nutrients \n\n',
         f'  0.0e-4   {dimensions[0] :.2e}  xlo xhi \n',f'  0.0e-4   {dimensions[1] :.2e}  ylo yhi \n',
         f'  0.0e-4   {dimensions[2] :.2e}  zlo zhi \n\n', ' Atoms \n\n'
         ]
    
    for i in range(1,n_cells+1):
        size = random.uniform(min_size, max_size)
        x = random.uniform(0,dimensions[0])
        y = random.uniform(0,dimensions[1])
        z = random.uniform(0,dimensions[2])
        L.append(f'     %d 1 {size :.2e}  1.3e3 {x :.2e} {y :.2e} {z :.2e} {size :.2e} \n'% (i))
    L.append('\n')
    L.append(' Nutrients \n\n')
    for i,nute in enumerate(Nutrients.keys()):
        if i == 0:
            L.append(f'     %d {nute} g {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} \n'% (i+1))
        else:
            L.append(f'     %d {nute} l {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} \n'% (i+1))
    L.append('\n')
    L.append(' Type Name \n\n')
    L.append('     1 ' + atomType + ' \n\n')
    L.append(' Diffusion Coeffs \n\n')
    for key in Diff_c.keys():
        L.append(f'     {key} {Diff_c[key]} \n')
    L.append('\n')
    L.append(' Growth Rate \n\n')
    L.append(f'     cyano {growthRate} \n\n')
    L.append(' Ks \n\n')    
    k = '     cyano'
    for key in K_s.keys():
        k = k + ' ' + str(K_s[key])
    k = k + f' \n\n'
    L.append(k) 
    for key in Params.keys():
        L.append(' ' + key + f' \n\n')
        L.append('     ' + atomType + ' ' + str(Params[key]) + ' \n\n')
    f= open(f"atom_{n}.in","w+")
    f.writelines(L)