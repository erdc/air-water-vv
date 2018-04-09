#!/usr/bin/env python
import os
#os.chdir('/home/travis/build/erdc/proteus/air-water-vv/2d/benchmarks/wavesloshing')
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
#import wavesloshing_so
import wavesloshing
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools

from proteus.defaults import (load_physics as load_p,
                              load_numerics as load_n,
                              load_system as load_so)

modulepath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../2d/benchmarks/wavesloshing')
petsc_options = os.path.join(os.path.dirname(os.path.abspath(__file__)),"../inputTemplates/petsc.options.asm")

class TestWaveSloshingTetgen(TestTools.AirWaterVVTest):

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
        FileList = ['wavesloshing.xmf',
                    'wavesloshing.h5']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass


            
    def test_run(self):
        from petsc4py import PETSc
        pList = []
        nList = []
        so = load_so('wavesloshing_so',modulepath)
        for (p,n) in so.pnList:
            pList.append(load_p(p,modulepath))
            nList.append(load_n(n,modulepath))
            if pList[-1].name == None:
                pList[-1].name = p
        #so = wavesloshing_so
        so.name = "wavesloshing"
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
        ns.calculateSolution('wavesloshing')
        assert(True)

        
#     def test_validate(self):
#         # Reading file
#         filename='pointGauge_levelset.csv'
#         with open (filename, 'rb') as csvfile: 
#             data=csv.reader(csvfile, delimiter=",")
#             a=[]
#             time=[]
#             probes=[]
#             nRows=0
#             for row in data:
#                 if nRows!=0:
#                     time.append(float(row[0]))
#                 if nRows==0:              
#                     for i in row:              
#                         if i!= '      time':   
#                             i=float(i[14:24])  
#                             probes.append(i)   
#                 row2=[]
#                 for j in row:
#                     if j!= '      time' and nRows>0.:
#                         j=float(j)
#                         row2.append(j)
#                 a.append(row2)
#                 nRows+=1
#             # Taking phi at the left boundary    
#             phi=[]
#             for k in range(1,nRows):
#                 phi.append(a[k][1]+0.05)  
#             # Validation of the results
#             # Phi a the left boundary at last time step
#             Phi_f_Cal = phi[-1] 
#             Phi_f_Ana = -wavesloshing.eta(0.0,time[len(time)-1])+0.05 
#             err = 100*abs(Phi_f_Ana-Phi_f_Cal)/Phi_f_Ana
# 	    assert(err<2.0) # Error < 2.0%

if __name__ == '__main__':
    pass
