#!/usr/bin/env python
import os
#os.chdir('/home/travis/build/erdc/proteus/air-water-vv/2d/oscillating_cylinder')
#os.chdir('../2d/oscillating_cylinder_new')
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
#import tank
#import tank_so
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools
#import AnalysisTools as at
import math
#import importlib
from proteus import defaults

from proteus.defaults import (load_physics as load_p,
                              load_numerics as load_n,
                              load_system as load_so)


modulepath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../2d/oscillating_cylinder')
petsc_options = os.path.join(os.path.dirname(os.path.abspath(__file__)),"../inputTemplates/petsc.options.asm") 


class TestOscillatingCylinderTetgen(TestTools.AirWaterVVTest):

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
            i = 0
            while i < len(all):
                if i < len(all)-1:
                    if all[i+1][0]!='-':
                        print "setting", all[i].strip(),all[i+1]
                        OptDB.setValue(all[i].strip('-'),all[i+1])
                        i = i+2
                    else:
                        print "setting", all[i].strip(), "True"
                        OptDB.setValue(all[i].strip('-'),True)
                        i = i+1
                else:
                    print "setting", all[i].strip(), "True"
                    OptDB.setValue(all[i].strip('-'),True)
                    i= i+1
        so.tnList=[0.0,0.001]+[0.001 + i*0.01 for i in range(1,int(round(0.03/0.01))+1)]            
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('tank')
        assert(True)


#    def test_validate(self):
#        probes = 'circle2D.csv'
#        model_Data = at.readProbeFile(probes)
#        model_time = model_Data[1]
#        model_data = model_Data[2]
#        model_fy = model_data[:,7]
#        
#        # removing the first second of simulation (noise)
#        pos_model_t_one_sec = np.where(model_time>1)[0][0] 
#        model_time = model_time[pos_model_t_one_sec:]
#        model_fy = model_fy[pos_model_t_one_sec:]
#        
#        rho = 998.2 #density
#        v = 0.8 #velocity
#        w = 1 #width
#        d = 0.25 #diameter
#        A = w*d
#        coeff = 0.5 * rho * A * V**2
#        average_fy = np.mean(model_fy)
#
#        cl_0 = (model_fy-average_fy)/coeff
#
#        cl_1 = np.zeros(len(model_time))
#        for n in range(len(t_model)-6):
#            cl_1[n+6] = np.mean(cl_0[n:n+12])
#
#
#
#        exp_T= 1.309
#        exp_A= 4.203
#        
#
#        gap = 0.85 # it can be changed
#
#        model_t_cor = model_time+gap #correctec with the gap
#
#        pos_model_t_40_s = np.where(model_t_cor>40)[0][0]
#        pos_model_t_45_s = np.where(model_t_cor>45)[0][0]
#
#        model_results = zeroCrossing(model_time[pos_model_t_40_s:pos_model_t_45_s],
#                cl_1[pos_model_t_40_s:pos_model_t_45_s])
#        model_T = model_results[0]
#        model_A = model_results[1]
#
#        err_T = 100 * abs(exp_T - model_T)/exp_T
#        err_A = 100 * abs(exp_A - model_A)/exp_A
#
#        assert(err_T<2)
#        assert(err_A<8)


if __name__=='__main__':
    pass
