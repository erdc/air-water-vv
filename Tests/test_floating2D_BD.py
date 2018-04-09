import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
import floating2D_BD_so
import os
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools
import AnalysisTools as at
import math

class TestFloatingCaissonBDTetgen(TestTools.AirWaterVVTest):

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
        FileList = ['floating2D_BD.xmf',
                    'floating2D_BD.h5']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass
            
    def test_run(self):
        from petsc4py import PETSc
        pList = []
        nList = []
        for (p,n) in floating2D_BD_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = floating2D_BD_so
        so.name = "floating2D_BD"
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
        ns.calculateSolution('floating2D_BD')
        assert(True)

        
#    def test_validate(self):
#        probes = 'caisson2D.csv'
#        datalist = at.readProbeFile(probes)
#        time = datalist[1]
#        data = datalist[2]
#        rz = data[:,5]
#        alpha = []
#        for i in range(0,len(rz)):
#            alpha.append(rz[i]*180/math.pi)
#        alpha = np.array(alpha)
#        it = np.where(time>2.5)[0][0]
#        period = at.zeroCrossing(time[:it],alpha[:it],up=False)[0]
#        period_ref = 0.93
#        err = 100*abs(period_ref-period)/abs(period_ref)
#        assert(err<7.0)


if __name__ == '__main__':
    pass
