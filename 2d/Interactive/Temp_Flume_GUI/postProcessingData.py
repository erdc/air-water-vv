import matplotlib.pyplot as plt
import numpy as np

plt.close('all')
scale=21
test=34
testName = 'S'+str(scale)+'T'+str(test)
with open('Results/' + testName + '/forceHistory_p.txt') as file:
    pData = file.readlines()

pXaxis = []
pYaxis = []
for strings in pData:
    l = strings.split()
    pXaxis.append(float(l[0]))
    pYaxis.append(float(l[1]))

with open('Results/' + testName + '/forceHistory_v.txt') as file:
    pData = file.readlines()

vXaxis = []
vYaxis = []
for strings in pData:
    l = strings.split()
    vXaxis.append(float(l[0]))
    vYaxis.append(float(l[1]))


t=[x*(60.0/len(pXaxis)) for x in range(len(pXaxis))]

f1, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(t,np.array(pXaxis)+np.array(vXaxis))
axarr[0].set_title('Pressure + Viscuous Force in X-axis')
axarr[0].set_ylabel('Force (N)')

axarr[1].plot(t,np.array(pYaxis)+np.array(vYaxis))
axarr[1].set_title('Pressure + Viscuous Force in Y-axis')
axarr[1].grid(True)
axarr[1].set_xlabel('time (s)')


#plt.setp([a.get_xticklabels() for a in axarr[0, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in axarr[:, 1]], visible=False)

# plt.xlim([0,15.0])
for ax in axarr.reshape(-1):
    ax.grid(True)


plt.show()

#### Spectral Analysis



