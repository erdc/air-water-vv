import numpy as np
from IPython.core.display import HTML

def Create_Images(name):
    
    from tables import  openFile
    import matplotlib.tri as mtri
    
    filename=name+'_p.h5'
    archive = openFile(filename,'r')
    with open('TimeList.txt') as file:
        tnList = file.readlines()
    tnList = [float(x.strip()) for x in tnList]
    nodes = archive.getNode("/nodesSpatial_Domain0")
    x=nodes[:,0]
    y=nodes[:,1]
    dim=[x.max(),y.max()]
    elements = archive.getNode("/elementsSpatial_Domain0")
    triang = mtri.Triangulation(x, y, elements)
    xg = np.linspace(0, dim[0], 20)
    yg = np.linspace(0, dim[1], 20)
    xi, yi = np.meshgrid(xg,yg)
    def get_N_RGBCol(N=5):
        HSV_tuples = [(x*1.0/N, 0.5, 0.5) for x in xrange(N)]
        return HSV_tuples
    #colors=get_N_RGBCol(2*len(domain.segments))
    import os   
    checkpath=os.path.exists('./PostProcessing') 
    if not checkpath:
        os.mkdir("PostProcessing")
    for it,t in enumerate(tnList):
        plt.clf()
        phi = archive.getNode("/phi_t"+`it`)
        vof = archive.getNode("/vof_t"+`it`)
        wvof = np.ones(vof.shape,'d')
        wvof -= vof
        u = archive.getNode("/u_t"+`it`)
        v = archive.getNode("/v_t"+`it`)
        f.clear()
        plt.xlabel(r'z[m]')
        plt.ylabel(r'x[m]')
        #plt.xlim(domain.x[0]-0.1*domain.L[0],domain.x[0]+domain.L[0]+0.1*domain.L[0])    
        #for si,s in enumerate(domain.segments):
        #    plt.plot([domain.vertices[s[0]][0],
        #                 domain.vertices[s[1]][0]],
        #                [domain.vertices[s[0]][1],
        #                 domain.vertices[s[1]][1]],
        #                color=colors[domain.segmentFlags[si]-1],
        #                linewidth=2,
        #                marker='o')
        plt.tricontourf(x,y,elements,wvof*np.sqrt(u[:]**2 + v[:]**2))
        plt.tricontour(x,y,elements,phi,[0], linewidth=4)
        u_interp_lin = mtri.LinearTriInterpolator(triang, u[:])
        v_interp_lin = mtri.LinearTriInterpolator(triang, v[:])
        u_lin = u_interp_lin(xi, yi)
        v_lin = v_interp_lin(xi, yi)
        plt.streamplot(xg, yg, u_lin, v_lin,color='k')
        plt.title('T=%2.3f' % (t,))
        plt.axis('equal')
        plt.xlim((0,dim[0]))
        f.tight_layout()
        f.canvas.draw()
        plt.savefig('PostProcessing/phi%4.4d.png' % (it,),dpi=200)

def Create_Video(name):
    data_uri_mp4 = open(name+".mp4", "rb").read().encode("base64").replace("\n", "")
    video_tag = """<video controls>
    <source type ="video/mp4" src="data:video/mp4;base64,{mp4}"/>
    Your browser does not support the video tag
    </video>""".format(mp4=data_uri_mp4)
    HTML(data=video_tag)