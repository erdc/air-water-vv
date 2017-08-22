showTriangles=True
showAir=False
def plot_current_results():
    global t_0, plt0
    """Makes a blocking call to retrieve remote data and displays the solution mesh
    as a contour plot.
    
    Parameters
    ----------
    in_place : bool
        By default it calls clear_output so that new plots replace old ones.  Set
        to False to allow keeping of all previous outputs.
    """
    import numpy as np
    import matplotlib.tri as mtri
    import math
    global nn, x, y, vof, triangles, t, phi, u, v, cfl, p, dt
    load_simulation_globals()
    Vmax=np.amax(np.sqrt(u[:]**2 + v[:]**2))
    print "u_max={0:.3f}, v_max={1:.3f}, Vmax={2:.3f}, cfl={3:.3f}, dt={4:.10f}".format(np.amax(u[:]),
                                                           np.amax(v[:]),
                                                           Vmax,
                                                           np.asarray(cfl).max(),
                                                           dt
                                                           )
    if math.isnan(np.amax(u[:])):
        return plt0
    f.clear()
    triang = mtri.Triangulation(x, y, triangles)

    dim=[x.max(),y.max()]
    xg = np.linspace(x.min(), dim[0], 20)
    yg = np.linspace(y.min(), dim[1], 20)
    xi, yi = np.meshgrid(xg,yg)
    
    wvof = np.ones(vof.shape,'d')
    if showAir:
        wvof = vof
    else:
        wvof -= vof
            
    plt.xlabel(r'z[m]')
    plt.ylabel(r'x[m]')
   
    #domain = plant.domain
    #def get_N_RGBCol(N=5):
    #    HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in xrange(N)]
    #    return HSV_tuples
    #colors=get_N_RGBCol(2*len(domain.segments))
    #for si,s in enumerate(domain.segments):
    #        plt.plot([domain.vertices[s[0]][0],
    #                     domain.vertices[s[1]][0]],
    #                     [domain.vertices[s[0]][1],
    #                     domain.vertices[s[1]][1]],
    #                     color=colors[domain.segmentFlags[si]-1],
    #                     linewidth=2,
    #                     marker='o')

    if showTriangles==True:
        plt.triplot(triang, linewidth=0.5)
    plt.tricontourf(x,y,triangles,np.sqrt(u[:]**2 + v[:]**2))
    clb=plt.colorbar(orientation='horizontal')
    #clb = plt.colorbar()
    #clb.set_label('Vmax', labelpad=-40, y=1.05, rotation=0)
    clb.set_label('Vmax')
    #plt.tricontourf(x,y,triangles,wvof)
    plt.tricontour(x,y,triangles,phi,[0], linewidth=4)
    u_interp_lin = mtri.LinearTriInterpolator(triang, u[:])
    v_interp_lin = mtri.LinearTriInterpolator(triang, v[:])
    u_lin = u_interp_lin(xi, yi)
    v_lin = v_interp_lin(xi, yi)
    plt.streamplot(xg, yg, u_lin, v_lin,color='k')
    plt.title('T=%3.5e' % (t,))
    plt.xlim((x.min(),dim[0]))
    t_0=t
    f.canvas.draw()
    plt0=plt
    return plt

'''Some plotting and monitoring utilities
'''
def load_simulation_globals():
    """Put some variables we need in engine namespace.

    These can then be retrieved by clients for inspection, visualization, etc.
    """
    global nn, x, y, vof, triangles, t, phi, u, v, cfl, p, dt
    model_vof = ns.modelList[1].levelModelList[-1]
    model_ls = ns.modelList[2].levelModelList[-1]
    # save solution and grid data for plotting purposes
    cfl = ns.modelList[0].levelModelList[-1].q[('cfl',0)]
    dt = ns.systemStepController.dt_system
    x = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,0]
    y = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:,1]
    triangles = ns.modelList[0].levelModelList[-1].mesh.elementNodesArray
    p = ns.modelList[0].levelModelList[-1].u[0].dof
    u = ns.modelList[0].levelModelList[-1].u[1].dof
    v = ns.modelList[0].levelModelList[-1].u[2].dof
    vof = ns.modelList[1].levelModelList[-1].u[0].dof_last
    phi = ns.modelList[2].levelModelList[-1].u[0].dof_last
    nn = len(x)
    #print "p={0}, u={1}, v={2}, triangles={3}, vof={4}, phi={5}".format(len(p),len(u),len(v),len(triangles),len(vof),len(phi))
    t=ns.systemStepController.t_system
    cfl*=dt
    
def simulation_alive():
    """Return True if the simulation thread is still running on any engine.
    """
    return simulation_thread.is_alive()

def monitor_simulation(refresh=5.0):
    global t_0
    """Monitor the simulation progress and call plotting routine.

    Supress KeyboardInterrupt exception if interrupted, ensure that the last 
    figure is always displayed and provide basic timing and simulation status.

    Parameters
    ----------
    refresh : float
      Refresh interval between calls to retrieve and plot data.  The default
      is 5s, adjust depending on the desired refresh rate, but be aware that 
      very short intervals will start having a significant impact.

    """
    import datetime as dt, time
    if not simulation_alive():
        f.clear()
        plt=plot_current_results()
        f.canvas.draw()
        print 'Simulation has already finished, no monitoring to do.'
        error= True
        return error
    t_0=0.
    t0 = dt.datetime.now()
    try:
        while simulation_alive():

            plt=plot_current_results()
            error= False
            tmon = dt.datetime.now() - t0
            t_sim=ns.systemStepController.t_system
            print 'Monitored for: %s. at t=%12.5e' % (tmon,t_sim)
            time.sleep(refresh) # so we don't hammer the server too fast
    except (KeyboardInterrupt):#, error.TimeoutError):
        msg = 'Monitoring interrupted, simulation is ongoing!'
    else:
        if ns.systemStepController.converged():
            msg =  "\x1b[31mStep Failed at t=%12.5e \x1b[0m" % (ns.systemStepController.t_system)
            error= True
        else:
            msg = 'Simulation completed!'

    tmon = dt.datetime.now() - t0
    print msg
    print 'Monitored for: %s.' % tmon
    return error
    
    