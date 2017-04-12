#!/usr/bin/env python
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
import crump_weir_so
import crump_weir as cw
import os
import numpy as np
import collections as cll
import csv
from proteus.test_utils import TestTools

class TestCrumpWeirTetgen(TestTools.AirWaterVVTest):

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
        FileList = ['crump_weir.xmf',
                    'crump_weir.h5']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)
            else:
                pass
            
    def test_run(self):
        from petsc4py import PETSc
        pList = []
        nList = []
        for (p,n) in crump_weir_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = crump_weir_so
        so.name = "crump_weir"
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
        ns.calculateSolution('crump_weir')
        assert(True)

        
    def test_validate(self):
        # Reading probes into the file
        file_vof='column_gauge.csv'
        file_u='u_over_crest.csv'
        
        def readProbeFile(filename):
            with open (filename, 'rb') as csvfile:
                data=np.loadtxt(csvfile, delimiter=",",skiprows=1)
                time=data[:,0]
                data = data[:,1:]
                csvfile.seek(0)
                header = csvfile.readline()
                header = header.replace("time","")
                header = header.replace("[","")
                header = header.replace("]","")
                header = header.replace(","," ")
                header = header.split()
                probeType = []
                probex = []
                probey = []
                probez = []        
                for ii in range(0,len(header),4):
                    probeType.append(header[ii])
                    probex.append(float(header[ii+1]))
                    probey.append(float(header[ii+2]))
                    probez.append(float(header[ii+3]))
                probeCoord = zip(np.array(probex),np.array(probey),np.array(probez))
                datalist = [probeType,probeCoord,time,data]
                return datalist
            
        data_vof = readProbeFile(file_vof)
        data_u = readProbeFile(file_u)      
        # Solves the problem if velocity and vof have not the same number of time steps
        if len(data_u[2])==len(data_vof[2]) or len(data_u[2])<len(data_vof[2]):
            time = data_u[2]
        else:
            time = data_vof[2]
        # Extracts the datas from the function readProbeFile   
        U = []
        VOF = []
        water_level = []
        for i in range(0,len(data_u[3])):
            U.append(data_u[3][i])
        for i in range(0,len(data_vof[3])):
            VOF.append(data_vof[3][i][16])
            water_level.append(cw.tank_dim[1]-data_vof[3][i][16])
        U = np.array(U)   
        VOF = np.array(VOF)
        water_level = np.array(water_level)
        u_coord = np.array(data_u[1])
        # Changes the order of the velocity
        tmp = []
        u = []
        for i in range(0,len(U[0])):
            for j in range(0,len(U)):
                tmp.append(U[j][i])
            u.append(tmp)
            tmp = []
        u = np.array(u)
        # Creates an average velocity of the water at each time step over the crest
        vel = []
        for i in range(0,len(time)):
            for j in range(0,len(u_coord)):
                if u_coord[j][1]<water_level[i]:
                    tmp.append(u[j][i])
            vel.append(np.mean(tmp))
            tmp = []
        vel = np.array(vel)
        # Calculates the time-average discharge over the crest
        section = []
        Q = []
        for i in range(0,len(time)):
            section.append(water_level[i]-0.5)
            Q.append(section[i]*vel[i])
        # Validation of the result
        Q_th = 2.0175 #Theoretical discharge between 20 s and 30 s
        T = time.tolist()
        T_20 = T.index(20.0)
        T_30 = T.index(30.0)
        T_20_to_30 = []
        for i in range(T_20,T_30):
            T_20_to_30.append(Q[i])
        Q_pr = np.mean(T_20_to_30) #Discharge between 20 s and 30 s obtained with PROTEUS
        err = 100*abs(Q_th-Q_pr)/Q_th
        assert(err<2.0)
        
if __name__ == '__main__':
    pass
