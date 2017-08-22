def plot_domain():
    import ipympl
    import matplotlib.tri as mtri
    from matplotlib import pyplot as  plt
    import numpy as np
    from proteus import Context
    f.clear()
    ct=Context.get()
    domain=ct.domain
    domain.L = ct.tank_dim
    domain.x = (0.,0.,0.)
    waterLine_x=ct.waterLine_x
    waterLine_z=ct.waterLine_z
    outflow_level=ct.outflow_level
    plt.xlabel(r'z[m]')
    plt.ylabel(r'x[m]')
    colors = ['b','g','r','c','m','y','k','w']
    #plt.xlim(domain.x[0]-0.1*domain.L[0],domain.x[0]+domain.L[0]+0.1*domain.L[0])    
    for si,s in enumerate(domain.segments):
        plt.plot([domain.vertices[s[0]][0],
                  domain.vertices[s[1]][0]],
                 [domain.vertices[s[0]][1],
                  domain.vertices[s[1]][1]],
                 color=colors[domain.segmentFlags[si]-1],
                 linewidth=2,
                 marker='o')
    plt.axis('equal')
    xg = np.linspace(0, domain.L[0], 20)
    yg = np.linspace(0, domain.L[1], 20)
    xi, yi = np.meshgrid(xg,yg)
    phi = np.zeros(xi.shape)
    for i in range(20):
        for  j in range(20):
            phi[i,j]  =ct.signedDistance([xg[j],yg[i]],waterLine_x,waterLine_z,outflow_level)
    plt.contourf(xg,yg,phi)
    plt.contour(xg,yg,phi,levels=[0],linewidths=[5], colors='k')
    plt.title('Bathmetry and Initial SDF')
    f.canvas.draw()

