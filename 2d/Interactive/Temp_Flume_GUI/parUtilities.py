import numpy as np




def mpi_order(seq):
    """Return elements of a sequence ordered by MPI rank.

    The input sequence is assumed to be ordered by engine ID."""
    ranks = view['rank']
    rank_indices = np.argsort(ranks)
    return [seq[x] for x in rank_indices]


showTriangles = True
showAir = False


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
    view.apply_sync(load_simulation_globals)
    x = np.concatenate(mpi_order(view['x']))
    y = np.concatenate(mpi_order(view['y']))
    u = np.concatenate(mpi_order(view['u']))
    v = np.concatenate(mpi_order(view['v']))
    vof = np.concatenate(mpi_order(view['vof']))
    shifts = np.cumsum([0] + mpi_order(view['nn'])[:-1])
    flat_triangles = np.concatenate([tri + shift for tri, shift in zip(mpi_order(view['triangles']), shifts)])

    t = mpi_order(view['t'])[0]
    phi = np.concatenate(mpi_order(view['phi']))
    cfl = mpi_order(view['cfl'])[0]
    p = np.concatenate(mpi_order(view['p']))
    dt = mpi_order(view['dt'])[0]

    Vmax = np.amax(np.sqrt(u[:] ** 2 + v[:] ** 2))
    print "u_max={0:.3f}, v_max={1:.3f}, Vmax={2:.3f}, cfl={3:.3f}, dt={4:.5f}".format(np.amax(u[:]),
                                                                                       np.amax(v[:]),
                                                                                       Vmax,
                                                                                       np.asarray(cfl).max() * dt,
                                                                                       dt
                                                                                       )


def load_simulation_globals():
    """Put some variables we need in engine namespace.

    These can then be retrieved by clients for inspection, visualization, etc.
    """
    global nn, x, y, vof, triangles, t, phi, u, v, cfl, p, dt
    model_vof = ns.modelList[1].levelModelList[-1]
    model_ls = ns.modelList[2].levelModelList[-1]
    # save solution and grid data for plotting purposes
    cfl = ns.modelList[0].levelModelList[-1].q[('cfl', 0)]
    dt = ns.systemStepController.dt_system
    x = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:, 0]
    y = ns.modelList[0].levelModelList[-1].mesh.nodeArray[:, 1]
    triangles = ns.modelList[0].levelModelList[-1].mesh.elementNodesArray
    p = ns.modelList[0].levelModelList[-1].u[0].dof
    u = ns.modelList[0].levelModelList[-1].u[1].dof
    v = ns.modelList[0].levelModelList[-1].u[2].dof
    vof = ns.modelList[1].levelModelList[-1].u[0].dof_last
    phi = ns.modelList[2].levelModelList[-1].u[0].dof_last
    nn = len(x)
    # print "p={0}, u={1}, v={2}, triangles={3}, vof={4}, phi={5}".format(len(p),len(u),len(v),len(triangles),len(vof),len(phi))
    t = ns.systemStepController.t_system
    cfl *= dt


def simulation_alive():
    """Return True if the simulation thread is still running on any engine.
    """
    return any(view.apply_sync(lambda: simulation_thread.is_alive()))


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
        plot_current_results()
        print 'Simulation has already finished, no monitoring to do.'
        error = True
        return error
    t_0 = 0.
    t0 = dt.datetime.now()
    try:
        while simulation_alive():
            plot_current_results()
            error = False
            tmon = dt.datetime.now() - t0
            t_sim = view.apply_sync(lambda: ns.systemStepController.t_system)[0]
            print 'Monitored for: %s. at t=%12.5e' % (tmon, t_sim)
            time.sleep(refresh)  # so we don't hammer the server too fast
    except (KeyboardInterrupt):  # , error.TimeoutError):
        msg = 'Monitoring interrupted, simulation is ongoing!'
    else:
        t_sim = view.apply_sync(lambda: ns.systemStepController.t_system)[0]
        T_des = view['T'][0]
        if not view.apply_sync(lambda: ns.systemStepController.converged()):
            msg = "\x1b[31mStep Failed at t=%12.5e \x1b[0m" % (t_sim)
            error = True
        else:
            msg = 'Simulation completed!'

    tmon = dt.datetime.now() - t0
    print msg
    print 'Monitored for: %s.' % tmon
    return error


