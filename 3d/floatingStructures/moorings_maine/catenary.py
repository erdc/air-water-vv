#http://mathhelpforum.com/calculus/96398-catenary-cable-different-heights.html
import numpy as np
import sys
from numpy import *
from math import *
nitercount = 0

def get_array(x):
    if np.isscalar(x):
        x = np.array([x])
    else:
        x = np.asarray(x)
    return x


def catenary_tension_elastic2(P1, P2, L, w, EA=None, a0=1., tol=1e-10, maxit=1000):
    """
    Calculates the solution for elastic catenary of unstretched length L between two points P1 and P2.
    
    Parameters
    ----------
    :param P1: coordinates of anchor (2D) (!) must be lower and on the left of fairlead
    :param P2: coordinates of point 2 (2D)  (!) must be higher and on the right of anchor
    :param L: length of cable
    :param w: submerged weight of cable
    :param EA: axial stiffness of cable
    :param a0: initial guess of solution to catenary equation (a=H/w)
    :param tol: tolerance for finding solution
    :param maxit: maximum of iterations
    """
    d = np.abs(P2[0]-P1[0])  # horizontal distance between points
    h = np.abs(P2[1]-P1[1])  # vertical distance between points
    minL = np.sqrt(h**2+d**2)+tol  # minimum line length

    L = get_array(L)  # unstretched line length
    w = get_array(w)  # submerged weight
    EA = get_array(EA)  # axial stiffness

    Lt = np.sum(L)  # total unstretched line length
    Ls = np.zeros(len(L))  # lifted line length
    Lsu = np.zeros(len(L))  # unstretched lifted line length
    e = np.zeros(len(L))  # stretching
    et = np.sum(e)  # total stretching
    Le = Lt+et  # stretched line length

    diff = tol+1
    niter = 0

    # cable straight to seabed: 
    # find tension and stretching
    for i in reversed(range(len(L))):
        Lsu[i] = L[i]
        for j in range(i, len(L)):
            e[i] = (w[i]*Lsu[i]/2.+np.sum(w[:i]*Lsu[:i]))*Lsu[i]/EA[i]
            Ls[i] = Lsu[i]+e[i]
        if np.sum(Ls) >= h:
            Lhi_low = 0
            Lhi_high = L[i]
            while diff > tol:
                Lsu[i] = (Lhi_low+Lhi_high)/2.
                for j in range(i, len(L)):
                    e[i] = (w[i]*Lsu[i]/2.+np.sum(w[:i]*Lsu[:i]))*Lsu[i]/EA[i]
                    Ls[i] = Lsu[i]+e[i]
                if np.sum(Ls) > h:
                    Lhi_high = Lsu[i]
                elif np.sum(Ls) < h:
                    Lhi_low = Lsu[i]
                diff = np.abs(np.sum(Ls)-h)
    # check if cable straight to seabed is solution
    if np.sum(L-Lsu+Ls)+tol >= h+d:
        # no horizontal tension
        state = 'long'
        a = 0
        Ls = h
        Le = Lt+et
        Tf = [0, np.sum(w*Lsu)]
        Ta = [0, 0]
        y_coords = lambda x: 1.
        xs = 0
        x0 = 0
        angle = 0
    else:
        # check if line is partly or fully lifted
        f = lambda a: a*(np.cosh(d/a)-1)-h
        a = bisection(f, 1., 1e15, tol=tol, maxit=maxit)
        Ls0 = a*sinh(d/a)  # maximum line length to be fully lifted
        # get actual line length assuming it is fully lifted (from parameter a)
        H = a*np.sum(w*L)/Lt
        Va = 0
        for i in range(len(e)):
            e[i] = np.sqrt(H**2+(Va+np.sum(w[:i]*Lsu[:i])+w[i]*Lsu[i]/2.)**2)*Lsu[i]/EA[i]
        Ls1 = Lt+np.sum(e)
        if Ls1 > Ls0:  # partly lifted
            a, e, Lsu = partly_lifted_elastic(d=d, h=h, L=L, w=w, EA=EA, maxit=maxit, tol=tol)
            state = 'part'
            angle = 0
            x0 = a*np.arccosh(1+h/a)
        elif Ls1 <= Ls0:  # fully lifted
            a, e, angle = fully_lifted_elastic(d=d, h=h, L=L, w=w, EA=EA, maxit=maxit, tol=tol)
            x0 = d
            if a is not np.nan:
                state = 'full'
                Lsu = L
            else:  # assume line is straight
                angle = np.arctan(h/d)
                H, e = straight(d=d, h=h, L=L, w=w, EA=EA, maxit=maxit, tol=tol)
                a = H/np.sum(w*L)
                state = 'straight'
    return a, e, x0, angle, state
    

def catenary_tension_elastic(P1, P2, L, w, EA=None, a0=1., tol=1e-10, maxit=1000):
    """
    Calculates the solution for elastic catenary of unstretched length L between two points P1 and P2.
    
    Parameters
    ----------
    :param P1: coordinates of anchor (2D) (!) must be lower and on the left of fairlead
    :param P2: coordinates of point 2 (2D)  (!) must be higher and on the right of anchor
    :param L: length of cable
    :param w: submerged weight of cable
    :param EA: axial stiffness of cable
    :param a0: initial guess of solution to catenary equation (a=H/w)
    :param tol: tolerance for finding solution
    :param maxit: maximum of iterations
    """
    d = np.abs(P2[0]-P1[0])  # horizontal distance between points
    h = np.abs(P2[1]-P1[1])  # vertical distance between points
    minL = np.sqrt(h**2+d**2)+tol  # minimum line length

    L = get_array(L)  # unstretched line length
    w = get_array(w)  # submerged weight
    EA = get_array(EA)  # axial stiffness

    Lt = np.sum(L)  # total unstretched line length
    Ls = np.zeros(len(L))  # lifted line length
    Lsu = np.zeros(len(L))  # unstretched lifted line length
    e = np.zeros(len(L))  # stretching
    et = np.sum(e)  # total stretching
    Le = Lt+et  # stretched line length

    diff = tol+1
    niter = 0

    # cable straight to seabed: 
    # find tension and stretching
    for i in reversed(range(len(L))):
        Lsu[i] = L[i]
        for j in range(i, len(L)):
            e[i] = (w[i]*Lsu[i]/2.+np.sum(w[:i]*Lsu[:i]))*Lsu[i]/EA[i]
            Ls[i] = Lsu[i]+e[i]
        if np.sum(Ls) >= h:
            Lhi_low = 0
            Lhi_high = L[i]
            while diff > tol:
                Lsu[i] = (Lhi_low+Lhi_high)/2.
                for j in range(i, len(L)):
                    e[i] = (w[i]*Lsu[i]/2.+np.sum(w[:i]*Lsu[:i]))*Lsu[i]/EA[i]
                    Ls[i] = Lsu[i]+e[i]
                if np.sum(Ls) > h:
                    Lhi_high = Lsu[i]
                elif np.sum(Ls) < h:
                    Lhi_low = Lsu[i]
                diff = np.abs(np.sum(Ls)-h)
    # check if cable straight to seabed is solution
    if np.sum(L-Lsu+Ls)+tol >= h+d:
        # no horizontal tension
        state = 'long'
        a = 0
        Ls = h
        Le = Lt+et
        Tf = [0, np.sum(w*Lsu)]
        Ta = [0, 0]
        y_coords = lambda x: 1.
        xs = 0
    else:
        # check if line is partly or fully lifted
        f = lambda a: a*(np.cosh(d/a)-1)-h
        a = bisection(f, 1., 1e15, tol=tol, maxit=maxit)
        Ls0 = a*sinh(d/a)  # maximum line length to be fully lifted
        # get actual line length assuming it is fully lifted (from parameter a)
        H = a*np.sum(w*L)/Lt
        Va = 0
        for i in range(len(e)):
            e[i] = np.sqrt(H**2+(Va+np.sum(w[:i]*Lsu[:i])+w[i]*Lsu[i]/2.)**2)*Lsu[i]/EA[i]
        Ls1 = Lt+np.sum(e)
        if Ls1 > Ls0:  # partly lifted
            a, e, Lsu = partly_lifted_elastic(d=d, h=h, L=L, w=w, EA=EA, maxit=maxit, tol=tol)
            state = 'part'
            angle = 0
        elif Ls1 <= Ls0:  # fully lifted
            a, e, angle = fully_lifted_elastic(d=d, h=h, L=L, w=w, EA=EA, maxit=maxit, tol=tol)
            if a is not np.nan:
                state = 'full'
                Lsu = L
            else:  # assume line is straight
                angle = np.arctan(h/d)
                H, e = straight(d=d, h=h, L=L, w=w, EA=EA, maxit=maxit, tol=tol)
                a = H/np.sum(w*L)
                state = 'straight'
        
        # postprocessing
        Ls = Lsu+e  # lifted line length
        Lst = np.sum(Ls)
        Lsut = np.sum(Lsu)
        w_av = np.sum(w*Lsu)/Lsut  # average submerged weight of cable
        H = a*w_av*(Lsut/Lst)  # horizontal tension
        Va = H*np.tan(angle)  # vertical tension at anchor
        Vf = Va+np.sum(w*Lsu)  # vertical tension at fairlead
        Ta = [H, Va]  # tension at anchor
        Tf = [H, Vf]  # tension at fairlead
        
        if state == 'part':
            x0 = a*np.arccosh(1+h/a)
            y_coords = lambda x: 0 if x < d-x0 else a*(np.cosh((x-(d-x0))/a)-1)
            dydx = lambda x: np.sinh((x)/a)
            angleP2 = np.arctan(dydx(x0))
            xs = [0]+[a*np.arcsinh(np.sum(Ls[:i+1])/a)+(d-x0) for i in range(len(Ls))] 
        elif state == 'full':
            Ls_tot = np.sum(L+e)
            xx = 0.5*(a*np.log((Ls_tot+h)/(Ls_tot-h))-d)
            xy = 0.5*(a*np.log((Ls_tot+h)/(Ls_tot-h))+d)
            k = xx-P1[0]
            m = P2[1]-a*np.cosh(xy/a)
            scorr = a*np.sinh(k/a)
            dydx = lambda x: np.sinh((x+k)/a)
            angleP2 = np.arctan(dydx(d))
            xs = [0]+[a*np.arcsinh((np.sum(Ls[:i+1])+scorr)/a)-k for i in range(len(Ls))]
            y_coords = lambda x: a*(np.cosh((x+k)/a))+m  # y coordinates of points along the
        elif state == 'straight':
            y_coords = lambda x: h/d*x
            xs = [0]+[h/d*np.sum(Ls[:i+1]) for i in range(len(Ls))]
    print(state)
    return Tf, Ta, Ls, e, a, y_coords, xs
    
def fully_lifted_elastic(d, h, L, w, EA, tol=1e-6, maxit=1000):
    Ls_tot = Le = 0
    Lt = np.sum(L)  # total length of cable
    w_av = np.sum(w*L/Lt) # average weight of cable
    e = np.zeros(len(L))  # stretching of cable segments
    
    t_high = h/d
    t_low = 0

    diff = 1.
    niter = 0
    a = 1.
    while diff > tol and niter < maxit:
        niter += 1
        if Le > Ls_tot:
            t_high = t
        elif Le < Ls_tot:
            t_low = t
        t = (t_low+t_high)/2.
        angle = np.arctan(t)
        # transcendental equation
        g = lambda a: a*(np.cosh(d/a+np.arcsinh(t))-np.cosh(np.arcsinh(t)))-h
        dg = lambda a: np.cosh(d/a+np.arcsinh(t))-d/a*np.sinh(d/a+np.arcsinh(t))
        a = bisection(f=g, int1=tol, int2=100000, tol=tol, maxit=maxit)
        #a = newton_raphson(f=g, df=dg, x0=a, tol=tol, maxit=maxit)
        #if a is np.nan:
        #    a = bisection(f=g, int1=1., int2=100000, tol=tol, maxit=maxit)
        # get new total Ls from solution a
        Ls_tot = np.sqrt((2*a*sinh(d/(2*a)))**2+h**2)
        # get new stretching from solution a
        Ta = a*w_av/np.cos(angle)*Lt/Ls_tot
        Ha = Ta*np.cos(angle)
        Va = Ta*np.sin(angle)
        for i in range(len(e)):
            e[i] = np.sqrt(Ha**2+(Va+(np.sum(w[:i]*L[:i])+w[i]*L[i]/2.))**2)*L[i]/EA[i]
        et = np.sum(e)
        Le = Lt+et  # store new Ls value as calculated with stretching
        x0 = d
        diff = np.abs(Le-Ls_tot)
    return a, e, angle


def fully_lifted_elastic_nofloor(d, h, L, w, EA, tol=1e-6, maxit=1000):
    Ls_tot = Le = 0
    Lt = np.sum(L)  # total length of cable
    w_av = np.sum(w*L/Lt) # average weight of cable
    e = np.zeros(len(L))  # stretching of cable segments
    
    t_high = h/d
    t_low = 0

    diff = 1.
    niter = 0
    a = 1.
    while diff > tol and niter < maxit:
        niter += 1
        if Le > Ls_tot:
            t_high = t
        elif Le < Ls_tot:
            t_low = t
        t = (t_low+t_high)/2.
        angle = np.arctan(t)
        # transcendental equation
        g = lambda a: 2.*a*np.sinh(d/(2.*a))-np.sqrt(Ls_tot**2.-h**2.)
        dg = lambda a: 2.*np.sinh(d/(2.*a))-d*np.cosh(d/(2.*a))/a
        a = bisection(f=g, int1=tol, int2=100000, tol=tol, maxit=maxit)
        #a = newton_raphson(f=g, df=dg, x0=a, tol=tol, maxit=maxit)
        #if a is np.nan:
        #    a = bisection(f=g, int1=1., int2=100000, tol=tol, maxit=maxit)
        # get new total Ls from solution a
        Ls_tot = np.sqrt((2*a*sinh(d/(2*a)))**2+h**2)
        # get new stretching from solution a
        Ta = a*w_av/np.cos(angle)*Lt/Ls_tot
        Ha = Ta*np.cos(angle)
        Va = Ta*np.sin(angle)
        for i in range(len(e)):
            e[i] = np.sqrt(Ha**2+(Va+(np.sum(w[:i]*L[:i])+w[i]*L[i]/2.))**2)*L[i]/EA[i]
        et = np.sum(e)
        Le = Lt+et  # store new Ls value as calculated with stretching
        x0 = d
        diff = np.abs(Le-Ls_tot)
    return a, e, angle

def fully_lifted_rigid(d, h, L, tol=1e-6, maxit=1000):
    L = np.sum(L)
    g = lambda a: 2.*a*np.sinh(d/(2.*a))-np.sqrt(L**2.-h**2.)
    dg = lambda a: 2.*np.sinh(d/(2.*a))-d*np.cosh(d/(2.*a))/a
    a0 = bisection(f=g, int1=1e-6, int2=1e15, tol=tol, maxit=maxit)
    a1 = newton_raphson(f=g, df=dg, x0=1., tol=tol, maxit=maxit)      
    if np.isnan(a1) or a1 < 0:
        a = a0
    else:
        a = a1
    return a

def full_catenary(d, h, L, tol=1e-6, maxit=1000):
    L = np.sum(L)
    g = lambda a: 2.*a*np.sinh(d/(2.*a))-np.sqrt(L**2.-h**2.)
    dg = lambda a: 2.*np.sinh(d/(2.*a))-d*np.cosh(d/(2.*a))/a
    # dg = lambda a: np.arccosh(d/a)+d*np.cosh(d/(2.*a))/a
    a = bisection(f=g, int1=1e-6, int2=1e15, tol=tol, maxit=maxit)
    # a = newton_raphson(f=g, df=dg, x0=1., tol=tol, maxit=maxit)      
    return a
    
def partly_lifted_elastic(d, h, L, w, EA, tol=1e-6, maxit=1000):
    diff = 1.
    niter = 0
    #g = lambda a: a*(np.cosh(d/2./a)-1)-h
    #a = bisection(f=g, int1=1., int2=100000, tol=tol, maxit=maxit)
    a = 1.
    e = np.zeros(len(L))
    Ls = np.zeros(len(L))
    Lt = np.sum(L)
    while diff > tol and niter < maxit:
        niter += 1
        Ls_tot = h*np.sqrt(1+2*a/h)  # lifted line length
        Ls_tot_check = 0
        end_line = False
        for i in reversed(range(len(L))):
            Ls_tot_check += L[i]+e[i]
            if Ls_tot_check < Ls_tot:
                Ls[i] = L[i]+e[i]
            elif Ls_tot_check > Ls_tot and end_line is False:
                Ls[i] = Ls_tot-np.sum(Ls[i+1:])  # find stretching under own weight
                end_line = True
            elif end_line is True:
                Ls[i] = 0
                e[i] = 0
        Lsu = Ls-e  # unstretched lifted line length
        Lsu_tot = np.sum(Lsu)
        w_av = np.sum(w*Lsu/Lsu_tot)
        x0 = a*np.arccosh(1+h/a)
        Ha = a*w_av*Lsu_tot/Ls_tot
        for i in range(len(e)):
            e[i] = np.sqrt(Ha**2+(np.sum(w[:i]*Lsu[:i])+w[i]*Lsu[i]/2.)**2)*Lsu[i]/EA[i]
        Le = Lt+np.sum(e)
        X0 = (Le-Ls_tot)+x0
        if Ls_tot > Le:
            diff = 1  # cannot find solution; line must be fully lifted
        else:
            diff = np.abs(X0-d)
        if diff > tol:
            a = a*((d/X0)**(d/x0))
    if niter >= maxit:
        a = e = Lsu = np.nan
    return a, e, Lsu

def partly_lifted_elastic2(d, h, L, w, EA, tol=1e-6, maxit=1000):
    diff = 1.
    niter = 0
    a = 1.
    e = np.zeros(len(L))
    Ls = np.zeros(len(L))
    Lt = np.sum(L)
    x0_high = d
    x0_low = 0
    Ls_tot = 0
    Lsu_tot = 0
    et = 0
    Lse = 0
    Ls = 0
    Lsu = np.zeros(len(L))
    while diff > tol and niter < maxit:
        niter += 1
        if Ls < Lse:
            x0_high = x0
        elif Ls > Lse:
            x0_low = x0
        x0 = (x0_low+x0_high)/2.
        lifted = False
        g = lambda a: a*(np.cosh(x0/a)-1.)-h
        dg = lambda a: np.cosh(d/a+np.arcsinh(t))-d/a*np.sinh(d/a+np.arcsinh(t))
        a = bisection(f=g, int1=1., int2=100000, tol=tol, maxit=maxit)
        Ls = h*np.sqrt(1+2*a/h)
        Lns_tot_check = 0
        ground = d-x0
        for i in (range(len(L))):
            if lifted is False:
                Lsu[i] = 0
                Lns_tot_check += L[i]
                if Lns_tot_check > ground:
                    Lsu[i] = Lns_tot_check-ground
                    lifted = True
            else:
                Lsu[i] = L[i]
        Lsu_tot = np.sum(Lsu)
        w_av = np.sum(w[i]*Lsu[i])/Ls
        H = a*w_av
        for i in range(len(L)):
            e[i] = np.sqrt(H**2+(np.sum(w[i:]*Lsu[i:])+w[i]*Lsu[i]/2.)**2)*Lsu[i]/EA[i]
        et = np.sum(e)
        Lse = Lsu_tot+et
        diff = np.abs(Ls-Lse)
    return a, e, Lsu
    

def partly_lifted_rigid(d, h, L, tol=1e-6, maxit=1000):
    diff = 1.
    niter = 0
    a = 1.
    Ls = np.zeros(len(L))
    Lt = np.sum(L)
    while diff > tol and niter < maxit:
        niter += 1
        Ls_tot = h*np.sqrt(1+2*a/h)  # lifted line length
        x0 = a*np.arccosh(1+h/a)
        X0 = (Lt-Ls_tot)+x0
        if Ls_tot > L:
            diff = 0  # cannot find solution; line must be fully lifted
            a = Ls = np.nan
        else:
            diff = np.abs(X0-d)
        if diff > tol:
            a = a*((d/X0)**(d/x0))    
    return a, Ls

def straight_elastic(d, h, L, w, EA, H_low=0, H_high=1e10, tol=1e-6, maxit=1000):

    Lt = np.sum(L)  # total length of cable
    w_av = np.sum(w*L/Lt) # average weight of cable
    e = np.zeros(len(L))  # stretching of cable segments
    et = (d**2+h**2)**0.5-Lt  # stretching to reach minimum length
    angle = np.arctan(h/d)  # angle
    
    H = (H_low+H_high)/2.  # guess of horizontal tension at anchor
    diff = 1.
    niter = 0
    while diff > tol and niter < maxit:
        niter += 1
        if et_check > et:
            H_high = H
        elif et_check < et:
            H_low = H
        H = (H_low+H_high)/2.
        Va = H*np.tan(angle)  # tension at anchor
        for i in range(len(e)):
            e[i] = ((H**2+(Va+w[i:]*L[w:]+w[i]*L[i]/2.)**2)**0.5)*L[i]/EA[i]
        et_check = np.sum(e)
    return H, e   


class MooringLine:
    """
    Class to create a mooring line between anchor and fairlead points
    :param w: submerged weight [N/m]
    :param E: Young's Modulus
    :param A0: cross-sectional area of unstretched cable [m^2]
    :param L: unstretched line length [m]
    :param UTS: ultimate tensile strength [N/mm^2]
    :param anchor: anchor coordinates
    :param fairlead: fairlead coordinates
    """
    count = 0
    def __init__(self, L, w, EA, UTS=None, anchor=None, fairlead=None, nd=3, tol=1e-6, maxiter=10000):
        self.__class__.count += 1
        self.nd = nd = nd
        self.name = 'mooring_line'+str(self.count)
        self.w = get_array(w)
        #self.E = E
        #self.A0 = A0
        self.EA = get_array(EA)
        self.L = get_array(L)
        self.UTS = UTS
        self.tol = tol
        self.maxiter = maxiter
        if anchor is None:
            self.anchor = np.zeros(nd)  # coordinates of anchor
        else:
            self.anchor = np.array(anchor)
        if fairlead is None:
            self.fairlead = np.zeros(nd)  # coordinates of fairlead
        else:
            self.fairlead = np.array(fairlead)
        self.anchor_coords_system = np.eye(3)
        self.setDirectionDistance()
        self.X = np.zeros(nd)  # local coords with anchor at origin
        self.x = np.zeros(nd)  # local coords (Ls, y)
        self.Tf = np.zeros(nd)  # horizontal pretension
        self.Ta = np.zeros(nd)
        self.Ls = 0.  # lifted line length
        self.e = 0.
        self.a = 1.
        self.x0 = None
        self.broken = False
        self.floor_level = 0.
        
        print(self.fairlead)
        self.setDirectionDistance()

    def setDirectionDistance(self):
        if self.nd == 3:
            self.distance = np.sqrt(np.sum((self.fairlead[:2]-self.anchor[:2])**2))
            self.direction = (self.fairlead[:2]-self.anchor[:2])/self.distance
        elif self.nd == 2:
            self.direction = np.array([1., 0.])
            self.distance = np.abs(self.fairlead[0]-self.anchor[0])
    def setVariables(self):
        if self.broken is False:
            if self.nd == 2:
                anchor2D = self.anchor
                fairlead2D = self.fairlead
                vec2D = np.array([1., 1.])
            elif self.nd == 3:
                #transf = np.linalg.inv(self.anchor_coords_system)
                #anchor_transf = np.dot(self.anchor, transf)
                #anchor2D = np.array([np.linalg.norm(anchor_transf[:2]), anchor_transf[2]])
                #fairlead_transf = np.dot(self.fairlead, transf)
                #fairlead2D = np.array([np.linalg.norm(fairlead_transf[:2]), fairlead_transf[2]])
                #vec = (fairlead_transf-anchor_transf)[:2]
                #vec = vec/np.linalg.norm(vec)
                anchor2D = np.array([0., 0.])
                self.setDirectionDistance()
                fairlead2D = np.array([self.distance, self.fairlead[2]-self.anchor[2]])
            if fairlead2D[0] < anchor2D[0]:
                # fairlead must be on the right of anchor
                fairlead2D[0], anchor2D[0] = anchor2D[0], fairlead2D[0]
                vec2D = np.array([-1., 1.])
        
  
        if self.floor_level is not None:
            a, e, x0, angle, state = catenary_tension_elastic2(P1=anchor2D, P2=fairlead2D, L=self.L, w=self.w, EA=self.EA, a0=self.a, tol=1e-8, maxit=10000)

        self.a = a
        self.e = e
        self.x0 = x0
        # postprocessing
        
        Ls = self.L+self.e
        Lsut = np.sum(self.L)
        Lst = np.sum(Ls)
        w_av = np.sum(self.w*self.L/Lsut)  # average submerged weight of cable
        H = a*w_av*(Lsut/Lst)  # horizontal tension
        print('H*w_av', a*w_av)
        Va = H*np.tan(angle)  # vertical tension at anchor
        Vf = Va+np.sum(self.w*self.L)  # vertical tension at fairlead
        Ta = [H, Va]  # tension at anchor
        Tf = [H, Vf]  # tension at fairlead

        h = np.abs(anchor2D[1]-fairlead2D[1])
        d = np.abs(anchor2D[0]-fairlead2D[0])
        Ls_tot = np.sum(self.L+self.e)
        angleP2 = 0
        if state == 'part':
            y_coords = lambda x: 0 if x < d-x0 else a*(np.cosh((x-(d-x0))/a)-1)
            dydx = lambda x: np.sinh((x)/a)
            # s_coords = lambda s: [s if s < d-x0 else a*np.arcsinh(s/a), 0 if s < d-x0 else a*(np.cosh((a*np.arcsinh(s/a)-(d-x0))/a)-1)]
            def s_coords(s):
                if s < d-x0:
                    return [s,0.]
                else:
                    s0 = (d-x0)
                    s = s-s0
                    x = s0+a*np.arcsinh(s/a)
                    y = a*(np.cosh((x-(d-x0))/a)-1)
                return [x, y]
            def ds_coords(s):
                if s < d-x0:
                    return [1.,0.]
                else:
                    s = s-(d-x0)
                    x = a/np.sqrt(a**2+s**2)
                    y = s/np.sqrt(a**2+s**2)
                return [x, y]
            angleP2 = np.arctan(dydx(x0))
            xs = [0]+[a*np.arcsinh(np.sum(Ls[:i+1])/a)+(d-x0) for i in range(len(Ls))] 
        elif state == 'full':
            xx = 0.5*(a*np.log((Ls_tot+h)/(Ls_tot-h))-d)
            xy = 0.5*(a*np.log((Ls_tot+h)/(Ls_tot-h))+d)
            k = xx-anchor2D[0]
            m = fairlead2D[1]-a*np.cosh(xy/a)
            scorr = a*np.sinh(k/a)
            dydx = lambda x: np.sinh((x+k)/a)
            angleP2 = np.arctan(dydx(d))
            xs = [0]+[a*np.arcsinh((np.sum(Ls[:i+1])+scorr)/a)-k for i in range(len(Ls))]
            y_coords = lambda x: a*(np.cosh((x+k)/a))+m  # y coordinates of points along the
            s_coords = lambda s: [a*np.arcsinh((scorr+s)/a)-k, a*np.cosh((a*np.arcsinh((scorr+s)/a)-k+k)/a)+m]
            def ds_coords(s):
                s = (s+scorr)
                x = a/np.sqrt(a**2+s**2)
                y = s/np.sqrt(a**2+s**2)
                return [x, y]
        elif state == 'straight':
            self.y_coords = lambda x: h/d*x
            #s_coords = lambda s: s if s < d-x0 else a*np.arcsinh(s/a), 0 if s < d-x0 else a*(np.cosh((x-(d-x0))/a)-1)
            xs = [0]+[h/d*np.sum(Ls[:i+1]) for i in range(len(Ls))]
        else:
            y_coords = lambda x: 1.


        print("H: ", state, H)
        print("H2: ", Vf*np.tan(angleP2))
        if self.nd == 2:
            self.Tf = Tf*vec2D
            self.Ta = Ta*vec2D
            self.y = y_coords
            self.s = s_coords
        if self.nd == 3:
            Tf3D = np.zeros(3)                
            Ta3D = np.zeros(3)
            # horizontal component from 2D to 3D on x-y
            Ta3D[:2] = Tf[0]*self.direction
            Tf3D[:2] = Tf[0]*self.direction
            self.y = lambda x: y_coords(x)*self.direction
            self.s = lambda s: self.anchor+[s_coords(s)[0]*self.direction[0], s_coords(s)[0]*self.direction[1], s_coords(s)[1]]
            self.ds = lambda s: [ds_coords(s)[0]*self.direction[0], ds_coords(s)[0]*self.direction[1], ds_coords(s)[1]]
            #Tf3D[:2] = Tf[0]*vec
            #Ta3D[:2] = Ta[0]*vec
            #self.y = lambda x, y: y_coords(np.sqrt(x**2+y**2))*vec
            # vertical component to z-axis
            Tf3D[2] = Tf[1]
            Ta3D[2] = Ta[1]
            # transform back with coord system to align with gravity
            #self.Tf = np.dot(self.anchor_coords_system, Tf3D)
            #self.Ta = np.dot(self.anchor_coords_system, Ta3D)
            self.Tf = Tf3D
            self.Ta = Ta3D

    def setAnchorCoords(self, coords):
        """
        Sets coordinates of anchor
        :param coords: coordinates of anchor
        """
        self.anchor[:] = np.array(coords)
        self.setDirectionDistance()

    def setFairleadCoords(self, coords):
        """
        Sets coordinates of fairlead
        :param coords: coordinates of fairlead
        """
        self.fairlead[:] = np.array(coords)
        self.setDirectionDistance()
        print('variables', self.distance, self.direction)
    
    def setCoords(self):
        transf = np.linalg.inv(self.anchor_coords_system)
        anchor = np.dot(self.anchor, transf)
        fairlead = np.dot(self.anchor, transf)
        
    def getTension(self, tol=1e-6, maxit=1000):
        """
        Gives the tension using line length and the coordinates of the anchor and fairlead
        :param tol: tolerance for calculation of tension
        :param maxit: maximum of iterations to calculate tension
        (!) works for seabed flat and perpendicular to gravity
        """
        if self.broken is False:
            if self.nd == 2:
                anchor2D = self.anchor
                fairlead2D = self.fairlead
                vec2D = np.array([1., 1.])
            elif self.nd == 3:
                transf = np.linalg.inv(self.anchor_coords_system)
                anchor_transf = np.dot(self.anchor, transf)
                anchor2D = np.array([np.linalg.norm(anchor_transf[:2]), anchor_transf[2]])
                fairlead_transf = np.dot(self.fairlead, transf)
                fairlead2D = np.array([np.linalg.norm(fairlead_transf[:2]), fairlead_transf[2]])
                vec = (fairlead_transf-anchor_transf)[:2]
                vec = vec/np.linalg.norm(vec)
            if fairlead2D[0] < anchor2D[0]:
                # fairlead must be on the right of anchor
                fairlead2D[0], anchor2D[0] = anchor2D[0], fairlead2D[0]
                vec2D = np.array([-1., 1.])

            Tf, Ta, self.Ls, self.e, self.a, y_func, xs = catenary_tension_elastic(P1=anchor2D, P2=fairlead2D, L=self.L, w=self.w, EA=self.EA, a0=self.a, tol=tol, maxit=maxit)

            if self.nd == 2:
                Tf = Tf*vec2D
                Ta = Ta*vec2D
                self.y = y_func
            if self.nd == 3:
                Tf3D = np.zeros(3)                
                Ta3D = np.zeros(3)
                # horizontal component from 2D to 3D on x-y
                Tf3D[:2] = Tf[0]*vec
                Ta3D[:2] = Ta[0]*vec
                self.y = lambda x, y: y_func(np.sqrt(x**2+y**2))*vec
                # vertical component to z-axis
                Tf3D[2] = Tf[1]
                Ta3D[2] = Ta[1]
                # transform back with coord system to align with gravity
                Tf = np.dot(self.anchor_coords_system, Tf3D)
                Ta = np.dot(self.anchor_coords_system, Ta3D)

            self.Tf = Tf
            self.Ta = Ta

            if self.UTS is not None:
                UBS = self.UTS*self.A0*(self.Ls-self.e)/self.Ls
                if np.linalg.norm(Tf) > UBS:
                    self.broken = True

def get_root_a(L, d, h, a0=1., tol=1e-6, maxit=1000):
    g = lambda a: 2.*a*np.sinh(d/(2.*a))-np.sqrt((L)**2.-h**2.)
    dg = lambda a: 2.*np.sinh(d/(2.*a))-d*np.cosh(d/(2.*a))/a
    a = newton_raphson(f=g, df=dg, x0=a0, tol=tol, maxit=maxit)      
    if np.isnan(a) or a < 0:
        a = bisection(f=g, int1=1., int2=1e15, tol=tol, maxit=maxit) 
    return a  
    
# Newton-Raphson
def newton_raphson(f, df, x0, tol=1e-6, maxit=1000):
    """
    Root finding algorithm (for transcendental equations)
    :param f: must be a function (so f(x) = 0 returns the required x)
    :param df: derivative of the function f (df/dx)
    :param x0: initial guess of x
    :param tol: tolerance
    :param maxit: maximum number of iterations
    :return: x
    """
    x_prev = x0
    x = x0-f(x0)/df(x0)
    err = np.abs(x-x_prev)
    niter = 0
    while err > tol and niter < maxit:
        niter += 1
        x_prev = x
        x = x-f(x)/df(x)
        err = np.abs(x-x_prev)
    #print('Newton-Raphson: iterations', niter, ', solution', x, ', err', err)
    if maxit <= niter:
        print('did not converge!')
        x = np.nan
    return x

def bisection(f, int1, int2, tol=1e-6, maxit=1000):
    """
    Root finding algorithm (for transcendental equations).
    A solution is found between a user-provided interval.
    :param f: must be a function (so f(x) = 0 returns the required x)
    :param int1: lower end value
    :param int2: higher end value
    :param tol: tolerance
    :param maxit: maximum number of iterations
    :return: x (solution)
    """
    err = np.abs(int2-int1)/2.
    niter = 0
    while err > tol and niter < maxit:
        niter += 1
        x = (int1+int2)/2.
        if np.sign(f(x)) == np.sign(f(int1)):
            int1 = x
        else:
            int2 = x
        err = np.abs(int2-int1)/2.
    #print('Bisection: iterations', niter, ', solution', x, ', err', err)
    if maxit <= niter:
        print('did not converge!')
        x = np.nan
    return x
