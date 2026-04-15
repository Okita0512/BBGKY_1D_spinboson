import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots() 
from matplotlib.pyplot import MultipleLocator, tick_params

lw = 2.0
legendsize = 32         # size for legend
font_legend = {'family':'Times New Roman', 'weight': 'roman', 'size':18}
color1 = 'black'            # 
color2 = 'r'

fs_to_au = 41.341

# ==============================================================================================
#                                      Fig 1a    
# ==============================================================================================

data = np.loadtxt("Pt.txt", dtype = float)                                  
plt.plot(data[:,0], data[:,1] - data[:,2], "--", linewidth = lw,  color = 'green', label = "tiers = 2")

plt.ylim(-1.0, 1.0)

plt.show()