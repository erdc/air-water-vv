#!/usr/bin/env python
import os
#os.chdir('/home/travis/build/erdc/proteus/air-water-vv/2d/numericalTanks/nonlinearWaves')
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
#import nonlinear_waves_so
#import nonlinear_waves as nlw
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools
#from AnalysisTools import readProbeFile,signalFilter,zeroCrossing,reflStat

from proteus.defaults import (load_physics as load_p,
                              load_numerics as load_n,
                              load_system as load_so)

modulepath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../2d/numericalTanks/nonlinearWaves')
petsc_options = os.path.join(os.path.dirname(os.path.abspath(__file__)),"../inputTemplates/petsc.options.asm")


class TestNonLinearWavesTetgen(TestTools.AirWaterVVTest):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        pass

    def teardown_method(self,method):
        """ Tear down function """
        FileList = ['nonlinear_waves.xmf',
                    'nonlinear_waves.h5']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass


            
    def test_run(self):
        from petsc4py import PETSc
        pList = []
        nList = []
        so = load_so('nonlinear_waves_so',modulepath)
        for (p,n) in so.pnList:
            pList.append(load_p(p,modulepath))
            nList.append(load_n(n,modulepath))
            if pList[-1].name == None:
                pList[-1].name = p
        so.name = "nonlinear_waves"
        if so.sList == []:
            for i in range(len(so.pnList)):
                s = default_s
                so.sList.append(s)
        Profiling.logLevel=7
        Profiling.verbose=True
        # PETSc solver configuration
        OptDB = PETSc.Options()
        with open(petsc_options) as f:
            all = f.read().split()
            i=0
            while i < len(all):
                if i < len(all)-1:
                    if all[i+1][0]!='-':
                        print "setting ", all[i].strip(), all[i+1]
                        OptDB.setValue(all[i].strip('-'),all[i+1])
                        i=i+2
                    else:
                        print "setting ", all[i].strip(), "True"
                        OptDB.setValue(all[i].strip('-'),True)
                        i=i+1
                else:
                    print "setting ", all[i].strip(), "True"
                    OptDB.setValue(all[i].strip('-'),True)
                    i=i+1
        so.tnList=[0.0,0.001,0.011]            
        #so.tnList=[0.0,0.001]+[0.001 + i*0.01 for i in range(1,int(round(0.03/0.01))+1)]            
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('nonlinear_waves')
        assert(True)

        
#     def test_validate(self):
#         # Reading probes into the file
#         file_vof = 'column_gauges.csv'
#         # Exctracting probes
#         data_vof = readProbeFile(file_vof)
#         time = data_vof[2]
#         vof = data_vof[3]
#         eta_num = []
#         tank_dim = nlw.opts.tank_dim
#         waterLevel = nlw.opts.water_level
#         i_mid = len(vof[0])/2-1
#         for i in range(0, len(vof)):
#             eta_num.append(tank_dim[1]-vof[i][i_mid]-waterLevel)
#         eta_num = np.array(eta_num)
#         # Theoretical eta
#         x = np.array(data_vof[1][2*i_mid])
#         wave = nlw.wave
#         eta_th = []
#         for i in range(0,len(time)):
#             eta_th.append(wave.eta(x,time[i]))
#         # Validation of the result
#         S = 0.
#         c = 0.
#         istart = np.where(time>=6.)[0][0]
#         iend = np.where(time>=18.)[0][0]
#         for i in range(istart,iend):
#             c = c + 1.
#             S = S + (eta_th[i]-eta_num[i])**2
#         err = np.sqrt(S/c)
#         err = 100*err/(nlw.opts.wave_height+waterLevel)
#         assert(err<5.0)
# 
#     def test_reflection(self):
#         dataW = readProbeFile('column_gauges.csv')
#         time = dataW[2]
#         L = nlw.opts.wave_wavelength
#         Nwaves = (nlw.opts.tank_dim[0]+nlw.opts.tank_sponge[0]+nlw.opts.tank_sponge[1])/L
#         T = nlw.opts.wave_period
#         Tend = time[-1]
#         Tstart = Tend-Nwaves*T
#         i_mid = len(dataW[3][0])/2-1
#         time_int = np.linspace(time[0],Tend,len(time))
#         data1 = np.zeros((len(time),len(dataW[3][0])),"d")
#         bf = 1.2
#         minf = 1./bf/T
#         maxf = bf / T
#         dx_array = nlw.opts.gauge_dx
#         Narray = int(round(L/6/dx_array))
#         data = np.zeros((len(data1),3))
#         zc = []
#         for ii in range(0,3):
#             data1[:,i_mid+ii*Narray] = np.interp(time_int,time,dataW[3][:,i_mid+ii*Narray])
#             data[:,ii] = signalFilter(time,data1[:,i_mid+ii*Narray],minf, maxf, 1.1*maxf, 0.9*minf)
#             zc.append(zeroCrossing(time,data[:,ii]))
#         H1 = zc[0][1]
#         H2 = zc[1][1]
#         H3 = zc[2][1]
#         HH = reflStat(H1,H2,H3,Narray*dx_array,L)[0]
#         RR = reflStat(H1,H2,H3,Narray*dx_array,L)[2]
#         assert(RR<0.3)
        
if __name__ == '__main__':
    pass
