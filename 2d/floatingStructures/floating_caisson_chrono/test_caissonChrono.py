import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
import os
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools
import AnalysisTools as at
import math
from proteus import Profiling

class TestFloatingCaissonChronoTetgen(TestTools.AirWaterVVTest):

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
        FileList = ['floating2D.xmf',
                    'floating2D.h5']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass

    fast = pytest.mark.skipif(not pytest.config.getoption("--runfast"), 
            reason="need --runfast option to run")
    
    slow = pytest.mark.skipif(pytest.config.getoption("--runfast"), 
            reason="no --runfast option to run")
    
    @fast
    def test_run_fast(self):
        os.chdir('2d/floatingStructures/floating_caisson_chrono')
        import floating2D_so
        import floating2D as f2d
        from petsc4py import PETSc
        pList = []
        nList = []
        for (p,n) in floating2D_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = floating2D_so
        so.name = "floating2D"
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
        so.tnList=[0.0,f2d.dt_init]+[f2d.dt_init + i*f2d.dt_out for i in range(1,int(round(0.3/f2d.dt_out)+1))]  
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('floating2D')
        assert(True)

    @slow
    def test_run_slow(self):
        import floating2D_so
        from petsc4py import PETSc
        pList = []
        nList = []
        for (p,n) in floating2D_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = floating2D_so
        so.name = "floating2D"
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
        ns.calculateSolution('floating2D')
        assert(True)
     
    @slow    
    def test_validate(self):
        probes = 'record_rectangle1.csv'
        datalist = at.readProbeFile(probes)
        time = datalist[1]
        data = datalist[2]
        rotq_e3 = data[:,7]
        alpha = []
        for i in range(0,len(rotq_e3)):
            alpha.append(2*math.asin(rotq_e3[i]))
            alpha[i] *= 360/(2*math.pi)
        alpha = np.array(alpha)
        it = np.where(time>2.5)[0][0]
        period = at.zeroCrossing(time[:it],alpha[:it],up=False)[0]
        period_ref = 0.93
        err = 100*abs(period_ref-period)/abs(period_ref)
        assert(err<4.0)


if __name__ == '__main__':
    pass
