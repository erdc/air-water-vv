#!/usr/bin/env python
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
import tank_so
import tank as tk
import os
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools

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
        FileList = ['tank.xmf',
                    'tank.h5']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass
            
    def test_run(self):
        from petsc4py import PETSc
        pList = []
        nList = []
        for (p,n) in tank_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = tank_so
        so.name = "tank"
        if so.sList == []:
            for i in range(len(so.pnList)):
                s = default_s
                so.sList.append(s)
        Profiling.logLevel=7
        Profiling.verbose=True
        # PETSc solver configuration
        OptDB = PETSc.Options()
        with open("../../../inputTemplates/petsc.options.asm") as f:
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
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('tank')
        assert(True)

        
    def test_validate(self):
        # Reading probes into the file
        file_vof = 'column_gauges.csv'

        def readProbeFile(filename):
            with open (filename, 'rb') as csvfile:
                data=np.loadtxt(csvfile, delimiter=",",skiprows=1)
                time=data[:,0]
                data = data[:,1:]
                csvfile.seek(0)
                header = csvfile.readline()
                header = header.replace("time","")
                header = header.replace("[","")
                header = header.replace("]","")
                header =  header.replace(","," ")
                header = header.split()
                probeType = []
                probex = []
                probey = []
                probez = []        
                for ii in range(0,len(header),4):
                    probeType.append(header[ii])
                    probex.append(float(header[ii+1]))
                    probey.append(float(header[ii+2]))
                    probez.append(float(header[ii+3]))
                probeCoord = zip(np.array(probex),np.array(probey),np.array(probez))
                datalist = [probeType,probeCoord,time,data]
                return datalist
            
        # Exctracting probes
        data_vof = readProbeFile(file_vof)
        time = data_vof[2]
        vof = data_vof[3]
        eta_num = []
        tank_dim = tk.opts.tank_dim
        waterLevel = tk.opts.water_level
        i_mid = len(vof[0])/2-1
        for i in range(0, len(vof)):
            eta_num.append(tank_dim[1]-vof[i][i_mid]-waterLevel)
        eta_num = np.array(eta_num)
        # Theoretical eta
        x = np.array(data_vof[1][2*i_mid])
        wave = tk.wave
        eta_th = []
        for i in range(0,len(time)):
            eta_th.append(wave.eta(x,time[i]))
        # Validation of the result
        S = 0.
        c = 0.
        istart = np.where(time>=6.)[0][0]
        iend = np.where(time>=18.)[0][0]
        for i in range(istart,iend):
            c = c + 1.
            S = S + (eta_th[i]-eta_num[i])**2
        err = np.sqrt(S/c)
        err = 100*err/(tk.opts.wave_height+waterLevel)
        assert(err<3.0)

        
        
if __name__ == '__main__':
    pass
