import pytest
import os
os.chdir('/home/travis/build/erdc/proteus/air-water-vv/2d/floatingStructures/floating_caisson_chrono')
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
#import floating2D_so
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools
import AnalysisTools as at
import math

from proteus.defaults import (load_physics as load_p,
                              load_numerics as load_n,
                              load_system as load_so)

modulepath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../2d/floatingStructures/floating_caisson_chrono')
petsc_options = os.path.join(os.path.dirname(os.path.abspath(__file__)),"../inputTemplates/petsc.option.asm")


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
            
    def test_run(self):
        from petsc4py import PETSc
        pList = []
        nList = []
        so = load_so('floating2D_so',modulepath)
        for (p,n) in so.pnList:
            pList.append(load_p(p,modulepath))
            nList.append(load_n(n,modulepath))
            if pList[-1].name == None:
                pList[-1].name = p
        #so = floating2D_so
        so.name = "floating2D"
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
        so.tnList=[0.0,0.001]+[0.001 + i*0.01 for i in range(1, int(round(0.03/0.01))+1)]            
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('floating2D')
        assert(True)

        
#    def test_validate(self):
#        probes = 'record_rectangle1.csv'
#        datalist = at.readProbeFile(probes)
#        time = datalist[1]
#        data = datalist[2]
#        rotq_e3 = data[:,7]
#        alpha = []
#        for i in range(0,len(rotq_e3)):
#            alpha.append(2*math.asin(rotq_e3[i]))
#            alpha[i] *= 360/(2*math.pi)
#        alpha = np.array(alpha)
#        it = np.where(time>2.5)[0][0]
#        period = at.zeroCrossing(time[:it],alpha[:it],up=False)[0]
#        period_ref = 0.93
#        err = 100*abs(period_ref-period)/abs(period_ref)
#        assert(err<4.0)


if __name__ == '__main__':
    pass
