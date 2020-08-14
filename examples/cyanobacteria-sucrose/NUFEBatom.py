import random
import argparse

parser = argparse.ArgumentParser(description='Create atom definition files')
parser.add_argument('--n', dest='num', action='store',
                   default=1,
                   help='Create atom definition files for NUFEB with --n #files desired (default is 1)')
parser.add_argument('--c',dest='culture_type',action='store',default='co',help='Set culture conditions with --c (default is co-culture)')
args = parser.parse_args()
if args.culture_type == 'co':
    cell_types = ['cyano','ecw']
    n_cyanos = int(random.uniform(1,100))
    n_ecw = int(random.uniform(1,100))
    n_cells = n_cyanos + n_ecw
elif args.culture_type == 'ax-c':
    cell_types = ['cyano']
    n_cyanos = int(random.uniform(1,100))
    n_cells = n_cyanos
elif args.culture_type == 'ax-e':
    cell_types = ['ecw']
    n_ecw = int(random.uniform(1,100))
    n_cells = n_ecw
growthRate = {'cyano' : 1.31e-5, 'ecw' : 2.7e-04} 
K_s = {'cyano' : {'sub' : 3.5e-4,'o2' : 2e-4, 'suc' : 1e-4,'co2' : 5e-1,'co2g' : 0},
       'ecw' : {'sub' : 0,'o2' : 1e-3, 'suc' : 3.4,'co2' : 5e-2,'co2g' : 0}}
Params = {'cyano' : {'Yield' : .55,'Maintenance' : 0,'Decay' : 0}, 
          'ecw' : {'Yield' : .43,'Maintenance' : 0,'Decay' : 0}}
for n in range(1,int(args.num)+1):
    dimensions = [1e-4,1e-4,1e-5]#x,y,z in meters
    Nutrients = {'sub' : 1e-1,'o2' : 9e-3, 'suc' : 1e-5, 'co2' : 4e-1, 'co2g' : 4e-1}
    NutesNum = len(Nutrients)
    Diff_c = {'sub' : 0,'o2' : 2.30e-9, 'suc' : 5.2e-10,'co2' : 1.9e-09, 'co2g' : 0}   
    
    L = [' NUFEB Simulation\r\n\n',f'     {n_cells} atoms \n',
         f'     {len(cell_types)} atom types \n',f'     {NutesNum} nutrients \n\n',
         f'  0.0e-4   {dimensions[0] :.2e}  xlo xhi \n',f'  0.0e-4   {dimensions[1] :.2e}  ylo yhi \n',
         f'  0.0e-4   {dimensions[2] :.2e}  zlo zhi \n\n', ' Atoms \n\n'
         ]
    if args.culture_type == 'co':
        atomType = 'cyano'
        min_size = 1.37e-6
        max_size = 1.94e-6
        for i in range(1,n_cyanos+1):
            size = random.uniform(min_size, max_size)
            x = random.uniform(0+size,dimensions[0]-size)
            y = random.uniform(0+size,dimensions[1]-size)
            z = random.uniform(0+size,dimensions[2]-size)
            L.append(f'     %d 1 {size :.2e}  375 {x :.2e} {y :.2e} {z :.2e} {size :.2e} \n'% (i))
        atomType = 'ecw'
        min_size = 8.8e-7
        max_size = 1.04e-6
        for i in range(n_cyanos+1,n_ecw+1+n_cyanos):
            size = random.uniform(min_size, max_size)
            x = random.uniform(0+size,dimensions[0]-size)
            y = random.uniform(0+size,dimensions[1]-size)
            z = random.uniform(0+size,dimensions[2]-size)
            L.append(f'     %d 2 {size :.2e}  375 {x :.2e} {y :.2e} {z :.2e} {size :.2e} \n'% (i))
    elif args.culture_type == 'ax-c':
        atomType = 'cyano'
        min_size = 1.37e-6
        max_size = 1.94e-6
        for i in range(1,n_cyanos+1):
            size = random.uniform(min_size, max_size)
            x = random.uniform(0+size,dimensions[0]-size)
            y = random.uniform(0+size,dimensions[1]-size)
            z = random.uniform(0+size,dimensions[2]-size)
            L.append(f'     %d 1 {size :.2e}  375 {x :.2e} {y :.2e} {z :.2e} {size :.2e} \n'% (i))
    elif args.culture_type == 'ax-e':
        atomType = 'ecw'
        min_size = 8.8e-7
        max_size = 1.04e-6
    
        for i in range(n_cyanos+1,n_ecw+1+n_cyanos):
            size = random.uniform(min_size, max_size)
            x = random.uniform(0+size,dimensions[0]-size)
            y = random.uniform(0+size,dimensions[1]-size)
            z = random.uniform(0+size,dimensions[2]-size)
            L.append(f'     %d 1 {size :.2e}  375 {x :.2e} {y :.2e} {z :.2e} {size :.2e} \n'% (i))
    L.append('\n')
    L.append(' Nutrients \n\n')
    for i,nute in enumerate(Nutrients.keys()):
        if 'sub' in nute:
            L.append(f'     %d {nute} g {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} \n'% (i+1))
        elif 'co2g' in nute:
            L.append(f'     %d {nute} g {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} \n'% (i+1))
        else:
            L.append(f'     %d {nute} l {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} {Nutrients[nute] :.2e} \n'% (i+1))
    L.append('\n')
    L.append(' Type Name \n\n')
    for i,cell in enumerate(cell_types,1):
        L.append(f'     {i} {cell} \n')
    L.append('\n')
    L.append(' Diffusion Coeffs \n\n')
    for key in Diff_c.keys():
        L.append(f'     {key} {Diff_c[key]} \n')
    L.append('\n')
    L.append(' Growth Rate \n\n')
    for cell in cell_types:
        L.append(f'     {cell} {growthRate[cell]} \n')
    L.append('\n')
    L.append(' Ks \n\n')   
    for cell in cell_types:
        k = f'     {cell}'
        for key in K_s[cell].keys():
            k = k + ' ' + str(K_s[cell][key])
        k = k + f' \n'
        L.append(k) 
    L.append('\n')

    for key in Params['cyano'].keys():
        L.append(' ' + key + f' \n\n')
        for cell in cell_types:
            L.append('     ' + cell + ' ' + str(Params[cell][key]) + ' \n')
            
    L.append('\n\n')
    f= open(f"atom_{n}.in","w+")
    f.writelines(L)