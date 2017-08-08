def plot_current_results():
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
    global nn, x, y, vof, elements, t, phi, u, v
    load_simulation_globals()
    
    triang = mtri.Triangulation(x, y, elements)
    
    domain = tank.domain
    xg = np.linspace(0, domain.L[0], 20)
    yg = np.linspace(0, domain.L[1], 20)
    xi, yi = np.meshgrid(xg,yg)
    
    wvof = np.ones(vof.shape,'d')
    wvof -= vof
    
    plt.xlabel(r'z[m]')
    plt.ylabel(r'x[m]')
    colors = ['b','g','r','c','m','y','k','w']
    for si,s in enumerate(domain.segments):
            plt.plot([domain.vertices[s[0]][0],
                         domain.vertices[s[1]][0]],
                        [domain.vertices[s[0]][1],
                         domain.vertices[s[1]][1]],
                        color=colors[domain.segmentFlags[si]-1],
                        linewidth=2,
                        marker='o')
            
    Vmax=np.amax(np.sqrt(u[:]**2 + v[:]**2))
    dt=tank.dt_fixed
    dx=tank.he
    cfl=Vmax*dt/dx
    if showTriangles:
        plt.triplot(triang, linewidth=0.5)
    plt.tricontourf(x,y,elements,wvof*(np.sqrt(u[:]**2 + v[:]**2)+0.2))
    plt.tricontour(x,y,elements,phi,[0], linewidth=4)
    u_interp_lin = mtri.LinearTriInterpolator(triang, u[:])
    v_interp_lin = mtri.LinearTriInterpolator(triang, v[:])
    u_lin = u_interp_lin(xi, yi)
    v_lin = v_interp_lin(xi, yi)
    plt.streamplot(xg, yg, u_lin, v_lin,color='k')
    
    plt.title('T=%2.2f, cfl=%2.2f, nodes=%d' % (t,cfl,len(elements)))
    plt.xlim((0,domain.L[0]))
    
    return plt

def load_simulation_globals():
    """Put some variables we need in engine namespace.

    These can then be retrieved by clients for inspection, visualization, etc.
    """
    global nn, x, y, vof, elements, t, phi, u, v
    model_vof = ns.modelList[1].levelModelList[-1]
    model_ls = ns.modelList[2].levelModelList[-1]
    model_twp = ns.modelList[0].levelModelList[-1]
    nodes = model_vof.mesh.nodeArray
    elements = model_vof.mesh.elementNodesArray
    x = nodes[:,0]
    y = nodes[:,1]
    vof = model_vof.u[0].dof_last
    phi = model_ls.u[0].dof_last
    u = model_twp.u[1].dof_last
    v = model_twp.u[2].dof_last
    nn = len(x)
    t=ns.systemStepController.t_system
    
def simulation_alive():
    """Return True if the simulation thread is still running on any engine.
    """
    return simulation_thread.is_alive()

def monitor_simulation(refresh=5.0):
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
        return
    
    t0 = dt.datetime.now()
    fig = None
    try:
        time.sleep(refresh) # so we don't hammer the server too fast
        while simulation_alive():
            f.clear()
            plt=plot_current_results()
            f.canvas.draw()
            tmon = dt.datetime.now() - t0
            print 'Monitored for: %s.' % tmon
            time.sleep(refresh) # so we don't hammer the server too fast
    except (KeyboardInterrupt):#, error.TimeoutError):
        msg = 'Monitoring interrupted, simulation is ongoing!'
    else:
        if ns.systemStepController.converged():
            msg =  "\x1b[31mStep Failed at t=%12.5e \x1b[0m" % (ns.systemStepController.t_system)
        else:
            msg = 'Simulation completed!'
    tmon = dt.datetime.now() - t0
    print msg
    print 'Monitored for: %s.' % tmon
    
    
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'