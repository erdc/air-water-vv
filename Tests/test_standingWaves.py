#!/usr/bin/env python
import os
#os.chdir('/home/travis/build/erdc/proteus/air-water-vv/2d/numericalTanks/standingWaves')
import pytest
from proteus.iproteus import *
#import standing_waves_so
#import standing_waves as sw
import numpy as np
import collections as cll
import csv
import math
from proteus.test_utils import TestTools

from proteus.defaults import (load_physics as load_p,
                              load_numerics as load_n,
                              load_system as load_so)

modulepath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../2d/numericalTanks/standingWaves')
petsc_options = os.path.join(os.path.dirname(os.path.abspath(__file__)),"../inputTemplates/petsc.options.asm")


class TestStandingWavesTetgen(TestTools.AirWaterVVTest):

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
        FileList = ['standing_waves.xmf',
                    'standing_waves.h5']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass


            
    def test_run(self):
        from petsc4py import PETSc
        pList = []
        nList = []
        so = load_so('standing_waves_so',modulepath)
        for (p,n) in so.pnList:
            pList.append(load_p(p,modulepath))
            nList.append(load_n(n,modulepath))
            if pList[-1].name == None:
                pList[-1].name = p
        #so = standing_waves_so
        so.name = "standing_waves"
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
        so.tnList=[0.0,0.001]+[0.001 + i*0.01 for i in range(1,int(round(0.03/0.01))+1)]            
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('standing_waves')
        assert(True)

        
#    def test_validate(self):
#        file_p = 'pressure_gaugeArray.csv'
#
#        def readProbeFile(filename):
#            with open (filename, 'rb') as csvfile:
#                data=np.loadtxt(csvfile, delimiter=",",skiprows=1)
#                time=data[:,0]
#                data = data[:,1:]
#                csvfile.seek(0)
#                header = csvfile.readline()
#                header = header.replace("time","")
#                header = header.replace("[","")
#                header = header.replace("]","")
#                header = header.replace(","," ")
#                header = header.split()
#                probeType = []
#                probex = []
#                probey = []
#                probez = []
#                for ii in range(0,len(header),4):
#                    probeType.append(header[ii])
#                    probex.append(float(header[ii+1]))
#                    probey.append(float(header[ii+2]))
#                    probez.append(float(header[ii+3]))
#                probeCoord = zip(np.array(probex),np.array(probey),np.array(probez))
#                datalist = [probeType,probeCoord,time,data]
#                return datalist
#
#        # Exctracting probes
#        data_p = readProbeFile(file_p)
#        T = sw.opts.wave_period #1.94
#        H = sw.opts.wave_height #0.025
#        depth = sw.opts.water_level #1.
#        L = sw.opts.wavelength #5.
#        time = data_p[2]
#        pressure = data_p[3][:,-1]
#        Z = -depth + data_p[1][0][1]
#        Nwaves = (sw.opts.tank_dim[0]+sw.opts.tank_sponge[0]+sw.opts.tank_sponge[1])/L
#
#        # Calculating the height with the pressure
#        def pressureToHeight(data,Z,depth,wavelength,rho,g):
#            k = 2*math.pi/wavelength
#            Kp = rho*g*np.cosh(k*(depth+Z))/np.cosh(k*depth)
#            return data/Kp
#        Tend = time[-1]
#        Tstart = Tend-Nwaves*T
#        Tend_period = Tstart+T
#        istart = np.where(time>=Tstart)[0][0]
#        iend_period = np.where(time>Tend_period)[0][0]
#        p_period = pressure[istart:iend_period]
#        p_0 = max(p_period)-min(p_period)
#        Hr = pressureToHeight(p_0,Z,depth,L,998.2,9.81)
#
#        # Validation of the result
#        err = 100*abs(2*H-Hr)/(2*H)
#        assert(err<12.0)
        
if __name__ == '__main__':
    pass
