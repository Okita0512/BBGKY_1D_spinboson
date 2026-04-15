import matplotlib.pyplot as plt
import numpy as np
fig, ax = plt.subplots() # 创建图实例
from matplotlib.pyplot import MultipleLocator, tick_params
from matplotlib import gridspec

lw = 2.0
legendsize = 32         # size for legend
font_legend = {'family':'Times New Roman', 'weight': 'roman', 'size':18}
color1 = 'black'            # 
color2 = 'r'


# ==============================================================================================
#                                      Fig 3a    
# ==============================================================================================

Unitlen = 6
lwid = 1.8
fig = plt.figure(figsize=(1.15 * Unitlen, 2.75 * Unitlen), dpi = 512)
spec = gridspec.GridSpec(ncols = 1, nrows = 3, height_ratios = [5, 5, 5])

ax0 = fig.add_subplot(spec[0])

data = np.loadtxt("Pt1_6.txt", dtype = float)
plt.plot(data[:,0], data[:,1] - data[:,2], ls = "--", color = 'green', label = 'tiers = 5', linewidth = lwid)

data = np.loadtxt("Pt1_11.txt", dtype = float)
plt.plot(data[:,0], data[:,1] - data[:,2], ls = "-", color = 'black', label = 'converged', linewidth = lwid)

data2 = np.loadtxt("DVR_1.txt", dtype = float)
plt.plot(data2[:,0], data2[:,1], ".", color = 'red', label = 'DVR')

# ==============================================================================================

time = 30.0             # x-axis range: (0, time)
y1, y2 = -1.2, 1.5     # y-axis range: (y1, y2)

# scale for major and minor locator
x_major_locator = MultipleLocator(5)
x_minor_locator = MultipleLocator(1)
y_major_locator = MultipleLocator(0.5)
y_minor_locator = MultipleLocator(0.1)

# x-axis and LHS y-axis
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.tick_params(which = 'major', length = 8, labelsize = 10)
ax.tick_params(which = 'minor', length = 4)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
y1_label = ax.get_yticklabels()
[y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params(labelsize = 20, which = 'both', direction = 'in')
plt.xlim(0.0, time)
plt.ylim(y1, y2)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 8)
ax2.tick_params(which = 'minor', length = 4)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(y1, y2)

# name of x, y axis and the panel
ax.set_xlabel('tΔ', font = 'Times New Roman', size = 18)
ax.set_ylabel('P(t)', font = 'Times New Roman', size = 18)
ax.legend(loc = 'upper center', prop = font_legend, markerscale = 1, frameon = False)
plt.legend(title = '(a)', frameon = False, title_fontsize = legendsize)

# ==============================================================================================
#                                      Fig 3b     
# ==============================================================================================

ax1 = fig.add_subplot(spec[1])

data = np.loadtxt("Pt2_11.txt", dtype = float)
plt.plot(data[:,0], data[:,1] - data[:,2], ls = "--", color = 'green', label = 'tiers = 10', linewidth = lwid)

# data = np.loadtxt("Pt2_16.txt", dtype = float)
# plt.plot(data[:,0], data[:,1] - data[:,2], ls = "--", color = 'blue', label = 'tiers = 15', linewidth = lwid)

data = np.loadtxt("Pt2_26.txt", dtype = float)
plt.plot(data[:,0], data[:,1] - data[:,2], ls = "-", color = 'black', label = 'converged', linewidth = lwid)

data2 = np.loadtxt("DVR_2.txt", dtype = float)
plt.plot(data2[:,0], data2[:,1], ".", color = 'red', label = 'DVR')

# ==============================================================================================

time = 20.0             # x-axis range: (0, time)
y1, y2 = -1.0, 1.8     # y-axis range: (y1, y2)

# scale for major and minor locator
x_major_locator = MultipleLocator(5)
x_minor_locator = MultipleLocator(1)
y_major_locator = MultipleLocator(0.5)
y_minor_locator = MultipleLocator(0.1)

# x-axis and LHS y-axis
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.tick_params(which = 'major', length = 8, labelsize = 10)
ax.tick_params(which = 'minor', length = 4)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
y1_label = ax.get_yticklabels()
[y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params(labelsize = 20, which = 'both', direction = 'in')
plt.xlim(0.0, time)
plt.ylim(y1, y2)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 8)
ax2.tick_params(which = 'minor', length = 4)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(y1, y2)

# name of x, y axis and the panel
ax.set_xlabel('tΔ', font = 'Times New Roman', size = 18)
ax.set_ylabel('P(t)', font = 'Times New Roman', size = 18)
ax.legend(loc = 'upper left', prop = font_legend, markerscale = 1, frameon = False)
plt.legend(title = '(b)', frameon = False, title_fontsize = legendsize)

# ==============================================================================================
#                                      Fig 3c 
# ==============================================================================================

ax2 = fig.add_subplot(spec[2])

data = np.loadtxt("Pt3_50.txt", dtype = float)
plt.plot(data[:,0], data[:,1] - data[:,2], ls = "--", color = 'green', label = 'tiers = 49', linewidth = lwid)

data = np.loadtxt("Pt3_75.txt", dtype = float)
plt.plot(data[:,0], data[:,1] - data[:,2], ls = "-", color = 'black', label = 'converged', linewidth = lwid)

data2 = np.loadtxt("DVR_3.txt", dtype = float)
plt.plot(data2[:,0], data2[:,1], ".", color = 'red', label = 'DVR')

# ==============================================================================================

time = 20.0             # x-axis range: (0, time)
y1, y2 = - 0.5, 1.0     # y-axis range: (y1, y2)

# scale for major and minor locator
x_major_locator = MultipleLocator(5)
x_minor_locator = MultipleLocator(1)
y_major_locator = MultipleLocator(0.5)
y_minor_locator = MultipleLocator(0.1)

# x-axis and LHS y-axis
ax = plt.gca()
ax.xaxis.set_major_locator(x_major_locator)
ax.xaxis.set_minor_locator(x_minor_locator)
ax.yaxis.set_major_locator(y_major_locator)
ax.yaxis.set_minor_locator(y_minor_locator)
ax.tick_params(which = 'major', length = 8, labelsize = 10)
ax.tick_params(which = 'minor', length = 4)

x1_label = ax.get_xticklabels()
[x1_label_temp.set_fontname('Times New Roman') for x1_label_temp in x1_label]
y1_label = ax.get_yticklabels()
[y1_label_temp.set_fontname('Times New Roman') for y1_label_temp in y1_label]

plt.tick_params(labelsize = 20, which = 'both', direction = 'in')
plt.xlim(0.0, time)
plt.ylim(y1, y2)

# RHS y-axis
ax2 = ax.twinx()
ax2.yaxis.set_major_locator(y_major_locator)
ax2.yaxis.set_minor_locator(y_minor_locator)
ax2.tick_params(which = 'major', length = 8)
ax2.tick_params(which = 'minor', length = 4)
ax2.axes.yaxis.set_ticklabels([])

plt.tick_params(which = 'both', direction = 'in')
plt.ylim(y1, y2)

# name of x, y axis and the panel
ax.set_xlabel('tΔ', font = 'Times New Roman', size = 18)
ax.set_ylabel('P(t)', font = 'Times New Roman', size = 18)

# legend location, font & markersize
ax.legend(loc = 'upper center', frameon = False, prop = font_legend) # lower right
# plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
plt.legend(title = '(c)', frameon = False, title_fontsize = legendsize)

plt.savefig("figure_1.eps", bbox_inches='tight')
# plt.show()
