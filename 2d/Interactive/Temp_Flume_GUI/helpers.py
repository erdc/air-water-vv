def plot_domain():
    import ipympl
    import matplotlib.tri as mtri
    from matplotlib import pyplot as  plt
    import numpy as np
    import Context 
    f.clear()
    ct=Context.get()
    domain=ct.domain
    domain.L = ct.tank_dim
    domain.x = (0.,0.,0.)
#    waterLine_x=ct.waterLine_x
    waterLine_z=ct.waterLine_z
#    outflow_level=ct.outflow_level
    plt.xlabel(r'z[m]')
    plt.ylabel(r'x[m]')
    import colorsys

    def get_N_RGBCol(N=5):
        HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in xrange(N)]
        rgb_out = []
        for rgb in HSV_tuples:
            rgb = map(lambda x: int(x*255),colorsys.hsv_to_rgb(*rgb))
            rgb_out.append(rgb)
        return HSV_tuples
    colors=get_N_RGBCol(2*len(domain.segments))

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
            phi[i,j]  =ct.signedDistance([xg[j],yg[i]],waterLine_z)
    plt.contourf(xg,yg,phi)
    plt.contour(xg,yg,phi,levels=[0],linewidths=[5], colors='k')
    plt.title('Bathmetry and Initial SDF')
    f.canvas.draw()

