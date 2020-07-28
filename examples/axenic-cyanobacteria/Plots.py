import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
SucRate = [float(s.split('_')[1]) for s in glob('./Sucrose*')]
types = ['./Sucrose_%.1f/Results/ntypes.csv'%i for i in SucRate]
biomass = ['./Sucrose_%.1f/Results/biomass.csv'%i for i in SucRate]
Cons = ['./Sucrose_%.1f/Results/avg_concentration.csv'%i for i in SucRate]
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
f, ax = plt.subplots()
for path,i in zip(types,SucRate):
    df = pd.read_csv(path,usecols=[0,1],names=['Time','Cells'],skiprows=1)
    ax.plot(df.Time*tStep2Days,df.Cells/CellNum2OD,label=i)
    # ax.plot(df.Time*tStep2Days,df.Cells*Biomass2OD,label=i)
ax.plot(Control.Time/24,Control.OD,label='Control')
ax.plot(IPTG.Time/24,IPTG.OD,label='Sucrose Export')
ax.set_xlabel('Time (days)')
ax.set_ylabel('OD')
ax.legend()

Suc = pd.DataFrame([],columns=['Time','Sucrose','Induction'])
for path,i in zip(Cons,SucRate):
    df = pd.read_csv(path,usecols=[0,2,3,4],names=['Time','O2','Sucrose','CO2(l)'],skiprows=1,index_col=0)
    df.index = df.index*tStep2Days
    df.O2 = df.O2/O2MW*1e3
    df.Sucrose = df.Sucrose/SucroseMW*1e3
    df['CO2(l)'] = df['CO2(l)']/CO2MW*1e3
    Suc = Suc.append(pd.DataFrame(pd.concat([pd.read_csv(path,usecols=[0,3],names=['Time','Sucrose'],skiprows=1),pd.Series(np.ones(len(df))*i,name='Induction')],axis=1),columns=['Time','Sucrose','Induction']),ignore_index=True)
    f, axes = plt.subplots()
    df.plot(ax=axes)
    axes.set_ylabel('Concentration (mM)')
    axes.set_xlabel('Time (days)')
    axes.set_title(f'Sucrose Ratio = {i}')
Suc.Sucrose = Suc.Sucrose/SucroseMW*1e3
