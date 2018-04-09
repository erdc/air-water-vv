#!/usr/bin/env python
import pytest
from proteus.iproteus import *
import submerged_breakwater_so
import submerged_breakwater as sbw
import os
import numpy as np
import collections as cll
import csv
import math
from proteus.test_utils import TestTools
from AnalysisTools import zeroCrossing

class TestSubmergedBreakwaterTetgen(TestTools.AirWaterVVTest):

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
        FileList = ['submerged_breakwater.xmf',
                    'submerged_breakwater.h5']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass
            
    def test_run(self):
        from petsc4py import PETSc
        pList = []
        nList = []
        for (p,n) in submerged_breakwater_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = submerged_breakwater_so
        so.name = "submerged_breakwater"
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
        ns.calculateSolution('submerged_breakwater')
        assert(True)

        
#    def test_validate(self):
#        file_vof = 'line_integral_gauges_1.csv'
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
#        data_vof = readProbeFile(file_vof)    
#        time = data_vof[2]
#        vof = data_vof[3]
#        ETA = []
#        tank_dim = sbw.tank_dim #[13.7815, 0.75] 
#        waterLevel = sbw.waterLevel #0.315
#        for j in range(len(data_vof[1])/2):
#            eta = []
#            for i in range(len(vof)):
#                eta.append(tank_dim[1]-vof[:,j][i]-waterLevel)
#            ETA.append(eta)
#        ETA = np.array(ETA)
#        zc = []
#        for i in range(len(ETA)):
#            zc.append(zeroCrossing(time,ETA[i]))
#        zc = np.array(zc)
#        K = np.mean(zc[2:][:,1])/sbw.opts.wave_height
#        Kref = 0.5
#        err = 100*abs(K-Kref)/Kref
#        assert(err<16.)
        
if __name__ == '__main__':
    pass
