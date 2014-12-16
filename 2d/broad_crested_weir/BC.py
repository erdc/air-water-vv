class boundaryConditions:

# Checking for entering the correct BC type
    def BCTypeCheck(self,BCType):
        boundaryList = ["pDirichlet","uDirichlet","vDirichlet", "wDirichlet","pAdvective","uAdvective","vAdvective","wAdvective","uDiffusive","vDiffusive","wDiffusive","vofDirichlet","vofAdvective"]
        if BCType not in boundaryList:
            print("Boundary condition type not valid")
            exit(1)


# Basic boundary conditions
    def empty(self):
        return lambda x,t: None    
    def constantValue(self,value):
        return lambda x,t: value



#Derived boundary conditions
    def noSlip(self,BCType):
        self.BCTypeCheck(BCType)
        if BCType == "uDirichlet" or BCType == "vDirichlet" or BCType == "wDirichlet" or ("Advective" in BCType):
            BC=self.constantValue(0.0)
            return BC
        else:
            BC=self.empty()
            return BC


                    
                
