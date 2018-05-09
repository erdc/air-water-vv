#!/usr/bin/env python
import os 
#os.chdir('/home/travis/build/erdc/proteus/air-water-vv/2d/caissonBreakwater/sliding')
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
#import tank_so
#import tank
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools
#import AnalysisTools as at
import math
from proteus import defaults

from proteus.defaults import (load_physics as load_p,
                              load_numerics as load_n,
                              load_system as load_so)


modulepath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../../2d/caissonBreakwater/sliding')
petsc_options = os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../inputTemplates/petsc.options.asm")


class NumericResults:

    @staticmethod
    def _read_log(file_name):
        log_file = open(file_name,'r')
        return log_file

    @classmethod
    def create_log(cls,file_name):
        sliding_log = cls._read_log(file_name)
        return sliding_log


class TestSlidingCaissonTetgen(TestTools.AirWaterVVTest):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        #pass
        Profiling.openLog("proteus.log",10)
        Profiling.logAllProcesses = True
            


    def teardown_method(self,method):
        Profiling.closeLog()

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
        so = load_so('tank_so',modulepath)
        for (p,n) in so.pnList:
            pList.append(load_p(p,modulepath))
            nList.append(load_n(n,modulepath))
            if pList[-1].name == None:
                pList[-1].name = p
        #so = tank_so
        so.name = "tank"
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
        so.tnList=[0.0,0.001,0.011]            
        #so.tnList=[0.0,0.001]+[0.001 + i*0.01 for i in range(1,int(round(0.03/0.01))+1)]            
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('tank')


        sliding_log = NumericResults.create_log('proteus.log')

        text = sliding_log.read()
        
        #if text.find('CFL') != -1:
        if text.find('Step Failed,') != -1:
   
            a = "No convergence"
        else:
            a = "good"
        

        if a  == "No convergence":
            print ("Convergence issue")
            assert False
        else:
            assert True

        
        
        #def failed(filename,word):
        #    file = open(filename,"r")
        #    text = file.read()

        #    if text.find(word) != -1:
        #        a = "No convergence"
        #    else:
        #        a = "good"
        #    file.close()
        #    return a

        #b = failed('proteus.log','Step Failed,')

        #if b == "No convergence":
        #    print ("Convergence issue")
        #    assert False
        #else:
        #    assert True

        #assert(True)


#    def test_validate(self):
#        probes = 'caisson2D.csv'
#        datalist = at.readProbeFile(probes)
#        time = datalist[1]
#        data = datalist[2]
#        ux_plastic = data[:,25]
#        t = []
#        j = 4 # time when the first wave gets the structure
#        wave_T = 1.3
#        ux_plast_per_wave = []
#        for n in range(13):
#            t.append(np.where(time>j)[0][0])
#            ux_plast_per_wave.append(ux_plastic[t[n]])
#            j = j + wave_T
#       
#        disp_ref = 0.0015
#        k = 0
#        disp = []
#        diff = []
#        for m in range (8): # 8 is the number of possible windows of 5 waves until the end of the simuation
#            disp.append(ux_plast_per_wave[k+5]-ux_plast_per_wave[k])
#            diff.append(abs(disp[m] - disp_ref))
#            k = k + 1
#        
#        pos_min = np.where(diff == min(diff))[0]
#
#        diff = np.array(diff)
#
#        err = 100 * (diff[pos_min]/disp_ref)
#        assert(err<10.0)
        

if __name__ == '__main__':
    pass





