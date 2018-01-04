import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
import tank_so
import os
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools
import AnalysisTools as at
import math

class TestSlidingCaissonTetgen(TestTools.AirWaterVVTest):

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
                        print "setting", all[i].strip(),all[i+1]
                        OptDB.setValue(all[i].strip('-'),all[i+1])
                        i=i+2
                    else:
                        print "setting", all[i].strip(), "True"
                        OptDB.setValue(all[i].strip('-'),True)
                        i=i+1
                else:
                    print "setting", all[i].strip(), "True"
                    OptDB.setValue(all[i].strip('-'),True)
                    i=i+1
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('tank')
        assert(True)


    def test_validate(self):
        probes = 'caisson2D.csv'
        datalist = at.readProbeFile(probes)
        time = datalist[1]
        data = datalist[2]
        ux_plastic = data[:,25]
        t = []
        j = 4 # time when the first wave gets the structure
        wave_T = 1.3
        ux_plast_per_wave = []
        for n in range(13):
            t.append(np.where(time>j)[0][0])
            ux_plast_per_wave.append(ux_plastic[t[n]])
            j = j + wave_T
       
        disp_ref = 0.0015
        k = 0
        disp = []
        diff = []
        for m in range (8): # 8 is the number of possible windows of 5 waves until the end of the simuation
            disp.append(ux_plast_per_wave[k+5]-ux_plast_per_wave[k])
            diff.append(abs(disp[m] - disp_ref))
            k = k + 1
        
        pos_min = np.where(diff == min(diff))[0]

        diff = np.array(diff)

        err = 100 * (diff[pos_min]/disp_ref)
        assert(err<10.0)
        

if __name__ == '__main__':
    pass





