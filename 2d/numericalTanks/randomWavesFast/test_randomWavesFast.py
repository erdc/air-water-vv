import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
import os
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools
from AnalysisTools import readProbeFile,signalFilter,zeroCrossing,reflStat

class TestRandomWavesFastTetgen(TestTools.AirWaterVVTest):

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
        FileList = ['random_waves.xmf',
                    'random_waves.h5']
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
        os.chdir('2d/numericalTanks/randomWavesFast')
        import random_waves_so
        import random_waves as rw
        from petsc4py import PETSc
        pList = []
        nList = []
        for (p,n) in random_waves_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = random_waves_so
        so.name = "random_waves"
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
        so.tnList=[0.0,rw.dt_init]+[rw.dt_init + i*rw.dt_out for i in range(1,int(round(0.1/rw.dt_out)+1))] 
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('random_waves')
        assert(True)

    @slow
    def test_run_slow(self):
        import random_waves_so
        from petsc4py import PETSc
        pList = []
        nList = []
        for (p,n) in random_waves_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = random_waves_so
        so.name = "random_waves"
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
        ns.calculateSolution('random_waves')
        assert(True)

    @slow    
    def test_validate(self):
        import random_waves as rw
        # Reading probes into the file
        file_vof = 'column_gauges.csv'
        # Exctracting probes
        data_vof = readProbeFile(file_vof)
        time = data_vof[2]
        vof = data_vof[3]
        # Eta for RandomWavesFast
        tank_dim = rw.opts.tank_dim
        waterLevel = rw.opts.water_level
        eta_fast = np.array(tank_dim[1]-vof[:,0]-waterLevel)
        # Eta for RandomWaves
        Tp = rw.Tp
        Hs = rw.Hs
        mwl = rw.mwl
        depth = rw.depth
        waveDir = np.array(rw.waveDir)
        g = np.array(rw.g)
        N = rw.N
        bandFactor = rw.bandFactor
        spectName = rw.spectName
        phi = rw.phi
        wave_ref = rw.wt.RandomWaves(Tp,Hs,mwl,depth,waveDir,g,N,bandFactor,spectName,spectral_params=None,phi=phi,fast=True)
        eta_ref = []
        X = np.array([0., 0., 0.])
        for i in range(0,len(time)):
            eta_ref.append(wave_ref.eta(X,time[i]))
        eta_ref = np.array(eta_ref)
        # Validation of the results
        S = 0.
        c = 0.
        for i in range(len(time)):
            c += 1.
            S += (eta_fast[i]-eta_ref[i])**2
        err = np.sqrt(S/c)
        err = 100*err/(rw.opts.Hs)
        assert(err<10.)


if __name__ == '__main__':
    pass
