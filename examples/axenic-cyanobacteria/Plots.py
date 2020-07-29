import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from glob import glob
SucRate = [float(s.split('_')[1]) for s in glob('./Sucrose*')]
types = ['./Sucrose_%.2f/Results/ntypes.csv'%i for i in SucRate]
biomass = ['./Sucrose_%.2f/Results/biomass.csv'%i for i in SucRate]
Cons = ['./Sucrose_%.2f/Results/avg_concentration.csv'%i for i in SucRate]
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
InitialCarbon = Volume*1.2e-1 #kg
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
for path,i,color in zip(types,SucRate,colors):
    df = pd.read_csv(path,usecols=[0,1],names=['Time','Cells'],skiprows=1)
    df.index = df.Time/60/60/24*10 #convert timesteps (10s) to days
    df.index.name='Days'
    df.iloc[:,1] = df.iloc[:,1]/CellNum2OD
    y0 = df.iloc[0,1]
    t = np.linspace(df.index[0], df.index[-1],1000)
    sol = odeint(monod_func, y0, t)
    ax.plot(df.iloc[:,1],label='NUFEB_suc' + str(i),c=color)
    ax.plot(t,sol,c=color,ls='--',label=r'Model_suc' + str(i))
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

#%%Nutrients
# Suc = pd.DataFrame([],columns=['Time','Sucrose','Induction'])
# for path,i in zip(Cons,SucRate):
#     df = pd.read_csv(path,usecols=[0,2,3,4],names=['Time','O2','Sucrose','CO2(l)'],skiprows=1,index_col=0)
#     df.index = df.index*tStep2Days
#     df.O2 = df.O2/O2MW*1e3
#     df.Sucrose = df.Sucrose/SucroseMW*1e3
#     df['CO2(l)'] = df['CO2(l)']/CO2MW*1e3
#     Suc = Suc.append(pd.DataFrame(pd.concat([pd.read_csv(path,usecols=[0,3],names=['Time','Sucrose'],skiprows=1),pd.Series(np.ones(len(df))*i,name='Induction')],axis=1),columns=['Time','Sucrose','Induction']),ignore_index=True)
#     f, axes = plt.subplots()
#     df.plot(ax=axes)
#     axes.set_ylabel('Concentration (mM)')
#     axes.set_xlabel('Time (days)')
#     axes.set_title(f'Sucrose Ratio = {i}')
# Suc.Sucrose = Suc.Sucrose/SucroseMW*1e3
