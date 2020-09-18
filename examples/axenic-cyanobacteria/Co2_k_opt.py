import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from glob import glob
import pickle
from sklearn.linear_model import LinearRegression
Run_folders = glob('./Run_*/')
def read_pkl(n):
# read python dict back from the file
    pkl_file = open(f"run_{n}.pkl", 'rb')
    f = pickle.load(pkl_file)
    pkl_file.close()
    return f
Run_params = [read_pkl(i) for i in range(1,len(Run_folders)+1)]

# K_s = []
# for i,run in enumerate(runs,1):
#     temp = open(f"atom_{i}.in","r").readlines()#
#     for j in range(len(temp)):
#         if '     cyano 0.00035 0.0002 0.01' in temp[j]:
#             # print(temp[j])
#             K_s.append(float(temp[j].split(' ')[-3]))
# K_s = [float(open(f"atom_{i}.in","r").readlines()[134].split(' ')[-2]) for i in range(1,20)]
types = ['./Run_%i/Results/ntypes.csv'%i for i in range(1,len(Run_folders)+1)]
biomass = ['./Run_%i/Results/biomass.csv'%i for i in range(1,len(Run_folders)+1)]
Cons = ['./Run_%i/Results/avg_concentration.csv'%i for i in range(1,len(Run_folders)+1)]
ExpPath = 'C:/Users/Jonathan/Documents/sucrose and growth CSCB-SPS.xlsx'
Control = pd.read_excel(ExpPath,sheet_name='Control')
IPTG = pd.read_excel(ExpPath,sheet_name='+IPTG')
SucroseMW = 342.3
O2MW = 32
CO2MW = 44.01
dens = 1e9
Volume = 1e-4*1e-4*1e-5 #m^3
CellNum2OD = Volume*1e6/0.3e-9
Biomass2OD = 1e12
tStep2Days = 10/3600/24
sucs = np.linspace(0,1,6)
pctC = 0.5
light = 1.00e-01 #kg/m^3
co2 = 4e-1
mu_max = 0.047*24 #2.25e-05# 1/d
K_m_light = 3.5e-04
K_m_co2 = 5e-2
# y= mu_max
colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a']
def monod_func(y,t):
    return y*mu_max * (light/(K_m_light + light)) * (co2/( K_m_co2 + co2))

# f.savefig('Growth_S9.png',dpi=600)
f, ax = plt.subplots(figsize=(14,9))
for path,i in zip(types,range(1,len(Run_folders)+1)):
    df = pd.read_csv(path,usecols=[0,1],names=['Time','Cells'],skiprows=1)
    df.index = df.Time/60/60/24*10 #convert timesteps (10s) to days
    df.index.name='Days'
    df.iloc[:,1] = df.iloc[:,1]/CellNum2OD
    y0 = df.iloc[0,1]
    t = np.linspace(df.index[0], df.index[-1],1000)
    sol = odeint(monod_func, y0, t)
    ax.plot(df.iloc[:,1],label=f'Sucrose Export {sucs[i-1] :.1f}')
    # ax.plot(t,sol,ls='--',label=f'monod {i:.2e}')
ax.set_yscale('log')
    # ax.set_xlabel('Time (D)')
    # ax.set_ylabel('Biomass (g)')
# ax.legend(frameon=False
ax.legend().remove()
    # ax.set_title('90% Sucrose Export')
    # ax.plot(df.Time*tStep2Days,df.Cells/CellNum2OD,label=i)
    # ax.plot(df.Time*tStep2Days,df.Cells*Biomass2OD,label=i)
    
ax.plot(Control.Time/24,Control.OD,label='Control',c='k')
ax.plot(IPTG.Time/24,IPTG.OD,label='Sucrose Export',c='grey')
ax.set_xlabel('Time (days)')
ax.set_ylabel('OD')
ax.legend()
#%%Carbon mass balance
viridis = plt.cm.get_cmap('viridis', 256)
color = viridis(np.linspace(0, 1, 20))
f, ax = plt.subplots(figsize=(14,9))
# plt.set_cmap(cmap=plt.get_cmap('viridis'))
j=0
for i,path in enumerate(biomass):
    InitialCarbon = Volume*12/CO2MW*Run_params[i]['Nutrients']['Concentration']['co2'] #kg
    df = pd.read_csv(path,usecols=[0,2],names=['Time','Biomass'],skiprows=1)
    df.index = df.Time/60/60/24*10 #convert timesteps (10s) to days
    df.index.name='Days'
    # df.iloc[:,1] = df.iloc[:,1]/(df.iloc[0,1]+InitialCarbon)
    # ax.plot(df.iloc[:,1]*pctC/(df.iloc[0,1]*pctC +InitialCarbon),label=f'K= {i:.2e}',c=color[j])
    ax.plot((df.iloc[:,1]*pctC - df.iloc[0,1]*pctC)/(InitialCarbon),label=f'Sucrose Export {sucs[i-1] :.1f}')
    j=j+1
    # ax.set_yscale('log')
ax.set_xlabel('Time (days)')
ax.set_ylabel('BiomassC/C')

ax.legend()

#%%Productivity
data = pd.DataFrame(columns=['SucroseRatio','Productivity'])

# Suc = pd.DataFrame([],columns=['Time','Sucrose','Induction'])
for cell_path,path,i in zip(types,Cons,range(1,len(Run_folders)+1)):
    df = pd.read_csv(cell_path,usecols=[0,1],names=['Time','Cells'],skiprows=1,index_col=0)
    df.index = df.index*tStep2Days*24
    df2 = pd.read_csv(path,usecols=[0,2,3,4],names=['Time','O2','Sucrose','CO2'],skiprows=1,index_col=0)
    df2.index = df2.index*tStep2Days*24

    df2.Sucrose = df2.Sucrose*Volume*1e18/df.Cells
    reg = LinearRegression()
    reg.fit(df2[:10].index.values.reshape(-1, 1),df2[:10].Sucrose)
    # Suc = Suc.append(pd.DataFrame(pd.concat([pd.read_csv(path,usecols=[0,3],names=['Time','Sucrose'],skiprows=1),pd.Series(np.ones(len(df))*i,name='Induction')],axis=1),columns=['Time','Sucrose','Induction']),ignore_index=True)
    f, axes = plt.subplots()
    df2.Sucrose.plot(ax=axes)
    axes.set_ylabel('Sucrose (fg)')
    axes.set_xlabel('Time (hrs)')
    axes.set_title(f'Sucrose Export {sucs[i-1] :.1f}')
    data = data.append(pd.DataFrame([[sucs[i-1],reg.coef_[0]]],columns=['SucroseRatio','Productivity']),ignore_index=True)
data.index = data.SucroseRatio
data.drop('SucroseRatio',axis=1,inplace=True)
f, ax = plt.subplots()
data.plot(ax=ax)
ax.set_ylabel('Sucrose Productivity (fg/cell/hr)')
ax.legend().remove()
#%% Nutrients
for path,i in zip(Cons,range(1,len(Run_folders)+1)):
    df = pd.read_csv(path,usecols=[0,2,3,4],names=['Time','O2','Sucrose','CO2'],skiprows=1,index_col=0)
    df.index = df.index*tStep2Days
    df.O2 = df.O2/O2MW*1e3
    df.Sucrose = df.Sucrose/SucroseMW*1e3
    df['CO2'] = df['CO2']/CO2MW*1e3
    # Suc = Suc.append(pd.DataFrame(pd.concat([pd.read_csv(path,usecols=[0,3],names=['Time','Sucrose'],skiprows=1),pd.Series(np.ones(len(df))*i,name='Induction')],axis=1),columns=['Time','Sucrose','Induction']),ignore_index=True)
    f, axes = plt.subplots()
    df.plot(ax=axes)
    axes.set_ylabel('Concentration (mM)')
    axes.set_xlabel('Time (days)')
    axes.set_title(f'Sucrose Export {sucs[i-1] :.1f}')
