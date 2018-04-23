#!/usr/bin/env python
import os
#os.chdir('/home/travis/build/erdc/proteus/air-water-vv/2d/benchmarks/dambreak_Ubbink')
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
#import dambreak_Ubbink_so
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools

from proteus.defaults import (load_physics as load_p,
                              load_numerics as load_n,
                              load_system as load_so)

modulepath = os.path.join(os.path.dirname(os.path.abspath(__file__)),'../2d/benchmarks/dambreak_Ubbink')
petsc_options = os.path.join(os.path.dirname(os.path.abspath(__file__)),"../inputTemplates/petsc.options.asm")


class NumericResults:

    @staticmethod
    def _read_log(file_name):
        log_file = open(file_name,'r')
        return log_file

    @classmethod
    def create_log(cls,file_name):
        ubbink_log = cls._read_log(file_name)
        return ubbink_log



class TestDambreakUbbinkTetgen(TestTools.AirWaterVVTest):

    @classmethod
    def setup_class(cls):
        pass

    @classmethod
    def teardown_class(cls):
        pass

    def setup_method(self,method):
        #pass
        Profiling.openLog("proteus.log",10)
        Profiling.logAllprocesses = True


    def teardown_method(self,method):
        Profiling.closeLog()

        """ Tear down function """
        FileList = ['dambreak_Ubbink.xmf',
                    'dambreak_Ubbink.h5']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass


            
    def test_run(self):
        from petsc4py import PETSc
        pList = []
        nList = []
        so = load_so('dambreak_Ubbink_so',modulepath)
        for (p,n) in so.pnList:
            pList.append(load_p(p,modulepath))
            nList.append(load_n(n,modulepath))
            if pList[-1].name == None:
                pList[-1].name = p
        #so = dambreak_Ubbink_so
        so.name = "dambreak_Ubbink"
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
        #so.tnList=[0.0,0.001]+[0.001 + i*0.01 for i in range(1,int(round(0.03/0.01))+1)]            
        ns = NumericalSolution.NS_base(so,pList,nList,so.sList,opts)
        ns.calculateSolution('dambreak_Ubbink')

        
        ubbink_log = NumericResults.create_log('proteus.log')

        text = ubbink_log.read()
        
        if text.find('dambreak_Colagrossi') != -1:
        #if text.find('Step Failed,') != -1:
   
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
        #
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
#        # Reading file
#        filename='pressureGauge.csv'
#        with open (filename, 'rb') as csvfile: 
#            data=csv.reader(csvfile, delimiter=",")
#            a=[]
#            time=[]
#            probes=[]
#            nRows=0
#            for row in data:
#                if nRows!=0:
#                    time.append(float(row[0]))
#                if nRows==0:              
#                    for i in row:              
#                        if i!= '      time':   
#                            i=float(i[14:24])  
#                            probes.append(i)   
#                row2=[]
#                for j in row:
#                    if j!= '      time' and nRows>0.:
#                        j=float(j)
#                        row2.append(j)
#                a.append(row2)
#                nRows+=1
#            # Takes the pressure   
#            pressure2=[]
#            for k in range(1,nRows):
#                pressure2.append(a[k][1])
#            # Validation of the results
#            maxPressure = max(pressure2)
#            s = 0
#            for i in range(1,len(pressure2)):
#                s = s+pressure2[i]
#            averagePressure = s/len(pressure2)
#            MaxPref = 3293.0
#            AvPref = 654.3
#            errMax = 100*abs(MaxPref-maxPressure)/MaxPref
#            errAv = 100*abs(AvPref-averagePressure)/AvPref
#            assert(errMax<2.0)
#            assert(errAv<0.5)

if __name__ == '__main__':
    pass
