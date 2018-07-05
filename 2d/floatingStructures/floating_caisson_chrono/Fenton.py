from __future__ import print_function
from __future__ import division
from past.utils import old_div
import numpy as np
from subprocess import check_call


def writeInput(waveheight, depth, period, mode='Period', current_criterion=1, current_magnitude=0, ncoeffs=8, height_steps=1, g=9.81):
    with open('Data.dat', 'w') as f:
        waveheight_dimless = old_div(waveheight,depth)
        mode = mode
        if mode == 'Period':
            length_dimless = period*np.sqrt(old_div(g,depth))
        elif mode == 'Wavelength':
            length_dimless = period*np.sqrt(old_div(g,depth))
        current_magnitude_dimless = old_div(current_magnitude,np.sqrt(g*depth))
        f.write('''Wave
{waveheight}
{mode}
{length}
{current_criterion}
{current_magnitude}
{ncoeffs}
{height_steps}
        '''.format(waveheight=waveheight_dimless, mode=mode, length=length_dimless, current_criterion=current_criterion,
        current_magnitude=current_magnitude_dimless, ncoeffs=ncoeffs, height_steps=height_steps))

def runFourier():
    check_call('./Fourier', shell=True)

def getBYCoeffs():
    FFT = np.genfromtxt('./Solution.res', delimiter=None, comments='#')
    BCoeffs = FFT[:,1]
    YCoeffs = FFT[:,2]
    return BCoeffs, YCoeffs

def copyFiles():
    from proteus import Profiling
    print(Profiling.logDir)
    import os
    check_call('cp ./Data.dat {0}'.format(os.path.join(Profiling.logDir, 'Data.dat')))
    check_call('cp ./Solution.res {0}'.format(os.path.join(Profiling.logDir, 'Solution.res')))
