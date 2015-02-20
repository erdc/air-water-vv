import numpy as np
import matplotlib as mpl
from matplotlib  import pyplot as plt
gauges = np.genfromtxt(fname="column_gauge.csv",dtype=('d','d','d'),delimiter=',',skip_header=1,names=("t","h1","h2"))
fig,(ax1,ax2) = plt.subplots(1,2)
ax1.plot(gauges["t"],gauges["h1"])
ax1.set_ybound(0.25,0.65)
ax2.plot(gauges["t"],gauges["h2"])
ax2.set_ybound(0.25,0.65)
plt.savefig("elevation_gauges.png")
