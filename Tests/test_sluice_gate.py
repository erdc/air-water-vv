#!/usr/bin/env python
import pytest
import os
#os.chdir('/home/travis/build/erdc/proteus/air-water-vv/2d/hydraulicStructures/sluice_gate')
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
#import sluice_gate_so
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools

from proteus.defaults import (load_physics as load_p,
                              load_numerics as load_n,
                              load_system as load_so)

modulepath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../2d/hydraulicStructures/sluice_gate')
petsc_options = os.path.join(os.path.dirname(os.path.abspath(__file__)),"../inputTemplates/petsc.options.asm")


class TestSluiceGateTetgen(TestTools.AirWaterVVTest):

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
        FileList = ['sluice_gate.xmf',
                    'sluice_gate.h5']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass
            
    def test_run(self):
        from petsc4py import PETSc
        pList = []
        nList = []
        so = load_so('sluice_gate_so',modulepath)
        for (p,n) in so.pnList:
            pList.append(load_p(p,modulepath))
            nList.append(load_n(n,modulepath))
            if pList[-1].name == None:
                pList[-1].name = p
        #so = sluice_gate_so
        so.name = "sluice_gate"
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
        #so.tnList=[0.0,0.001]+[0.001 + i*0.01 for i in range(1, int(round(0.03/0.01))+1)]            
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('sluice_gate')
        assert(True)

        
#    def test_validate(self):
#        # Reading probes into the file
#        filename='combined_column_gauge.csv'
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
#                header =  header.replace(","," ")
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
#        # Extracts the datas from the function readProbeFile  
#        datalist = readProbeFile(filename)
#        time = datalist[2]
#        # Calculates the time-average discharge under the gate
#        U = []
#        for i in range(0,len(datalist[3])):
#            U.append(np.mean(datalist[3][i]))
#        U = np.array(U)
#        Q = U*0.25
#        Q_th = 1.037 #Theoretical discharge between 20 s and 30 s
#        T = time.tolist()
#        T_20 = T.index(20.0)
#        T_30 = T.index(30.0)
#        T_20_to_30 = []
#        for i in range(T_20,T_30):
#            T_20_to_30.append(Q[i])
#        Q_pr = np.mean(T_20_to_30) #Discharge between 20 s and 30 s obtained with PROTEUS
#        err = 100*abs(Q_th-Q_pr)/Q_th
#        assert(err<7.)
        
if __name__ == '__main__':
    pass
