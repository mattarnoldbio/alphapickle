#%%
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm as cm
import numpy as np
import matplotlib.colors as col
from scipy.interpolate import interp1d

pLDDT = pd.read_csv("/Users/mattarnold/OneDrive - University of Glasgow/AlphafoldRequestResults/20212022_Alain/result_model_3_ptm_pLDDT.csv")
#%%
x = list(pLDDT['Unnamed: 0'])
y = list(pLDDT["lDDT"])
cmap = col.LinearSegmentedColormap.from_list("", ["red","orange","yellow","cornflowerblue","blue"])

plt.figure(figsize=(12,6))
ticks = np.arange(0, pLDDT['Unnamed: 0'].size, 100)
plt.xticks(ticks, fontname="Times New Roman")
plt.yticks(fontname="Times New Roman")
plt.xlabel("Residue index", size=14, fontweight="bold", fontname="Times New Roman")
plt.ylabel("Predicted LDDT",size=12, fontweight="bold", fontname="Times New Roman")
#plt.plot(x,y)
plt.scatter(x,y, c=y, cmap=cmap, s=5) #,edgecolor='none')
plt.clim(0,100)
scale = plt.colorbar(shrink=0.5)
scale.set_label(label="Predicted LDDT",size=12, fontweight="bold", fontname="Times New Roman")
plt.show

# %%
