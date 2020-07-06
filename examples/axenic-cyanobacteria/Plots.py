import pandas as pd
import matplotlib.pyplot as plt
SucRate = [0,0.1,0.3,0.5,0.7]
paths = ['./Sucrose_%.1f/Results/ntypes.csv'%i for i in SucRate]
f, ax = plt.subplots()
for path in paths:
    df = pd.read_csv(path,usecols=[0,1],names=['Time','Cells'],skiprows=1)
    ax.plot(df.Time,df.Cells)