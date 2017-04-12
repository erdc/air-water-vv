#!/usr/bin/env python
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
import wavesloshing_so
import wavesloshing
import os
import numpy as np
import collections as cll
import csv


class TestWaveSloshingTetgen():

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
        for (p,n) in wavesloshing_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = wavesloshing_so
        so.name = "wavesloshing"
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
        ns.calculateSolution('wavesloshing')
        assert(True)

        
    def test_validate(self):
        # Reading file
        filename='pointGauge_levelset.csv'
        with open (filename, 'rb') as csvfile: 
            data=csv.reader(csvfile, delimiter=",")
            a=[]
            time=[]
            probes=[]
            nRows=0
            for row in data:
                if nRows!=0:
                    time.append(float(row[0]))
                if nRows==0:              
                    for i in row:              
                        if i!= '      time':   
                            i=float(i[14:24])  
                            probes.append(i)   
                row2=[]
                for j in row:
                    if j!= '      time' and nRows>0.:
                        j=float(j)
                        row2.append(j)
                a.append(row2)
                nRows+=1
            # Taking phi at the left boundary    
            phi=[]
            for k in range(1,nRows):
                phi.append(a[k][1]+0.05)  
            # Validation of the results
            # Phi a the left boundary at last time step
            Phi_f_Cal = phi[-1] 
            Phi_f_Ana = -wavesloshing.eta(0.0,time[len(time)-1])+0.05 
            err = 100*abs(Phi_f_Ana-Phi_f_Cal)/Phi_f_Ana
	    assert(err<2.0) # Error < 2.0%

if __name__ == '__main__':
    pass
