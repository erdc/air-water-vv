#!/usr/bin/env python
import os
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools

class TestSharpCrestedWeirTetgen(TestTools.AirWaterVVTest):

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
        FileList = ['sharp_crested_weir.xmf',
                    'sharp_crested_weir.h5']
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
        os.chdir('2d/hydraulicStructures/sharp_crested_weir')
        import sharp_crested_weir_so
        from petsc4py import PETSc
        pList = []
        nList = []
        for (p,n) in sharp_crested_weir_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = sharp_crested_weir_so
        so.name = "sharp_crested_weir"
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
        so.tnList=[0.0,0.001]+[0.001 + i*0.02 for i in range(1,int(round(0.05/0.02)+1))]  
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('sharp_crested_weir')
        assert(True)
   
    @slow        
    def test_run_slow(self):
        import sharp_crested_weir_so
        from petsc4py import PETSc
        pList = []
        nList = []
        for (p,n) in sharp_crested_weir_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = sharp_crested_weir_so
        so.name = "sharp_crested_weir"
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
        ns.calculateSolution('sharp_crested_weir')
        assert(True)
        
       
if __name__ == '__main__':
    pass
