#!/usr/bin/env python
import pytest
from proteus.iproteus import *
from proteus import Comm
comm = Comm.get()
import dambreak_Colagrossi_so

class TestDambreakCollagrossiTetgen():

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
        FileList = ['dambreak_Colagrossi.xmf',
                    'dambreak_Colagrossi.h5']
        for file in FileList:
            if os.path.isfile(file):
                os.remove(file)

    def test_run(self):
        from petsc4py import PETSc
        pList = []
        nList = []
        for (p,n) in dambreak_Colagrossi_so.pnList:
            pList.append(__import__(p))
            nList.append(__import__(n))
            if pList[-1].name == None:
                pList[-1].name = p
        so = dambreak_Colagrossi_so
        so.name = "dambreak_Colagrossi"
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
        ns.calculateSolution('dambreak_Colagrossi')
        assert(True)

if __name__ == '__main__':
    pass
