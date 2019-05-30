def geometry_to_gmsh(domain):
    import py2gmsh
    from py2gmsh.Mesh import *
    from py2gmsh.Entity import *
    from py2gmsh.Field import *
    self = domain
    lines_dict = {}

    mesh = Mesh()

    if self.boundaryTags:
        for tag, flag in self.boundaryTags.items():
            phys = PhysicalGroup(nb=flag, name=tag)
            mesh.addGroup(phys)

    for i, v in enumerate(self.vertices):
        if domain.nd == 2:
            p = Point([v[0], v[1], 0.])
        else:
            p = Point(v)
        mesh.addEntity(p)
        g = mesh.groups.get(self.vertexFlags[i])
        if g:
            g.addEntity(p)
    nb_points = i+1
    for i in range(nb_points):
        lines_dict[i] = {}
        
    for i, s in enumerate(self.segments):
        lines_dict[s[0]][s[1]] = i
        l = Line([mesh.points[s[0]+1], mesh.points[s[1]+1]])
        mesh.addEntity(l)
        g = mesh.groups.get(self.segmentFlags[i])
        if g:
            g.addEntity(l)

    for i, f in enumerate(self.facets):
        if self.nd == 3 or (self.nd == 2 and i not in self.holes_ind):
            lineloops = []
            for j, subf in enumerate(f):
                lineloop = []
                # vertices in facet
                for k, ver in enumerate(subf):
                    if ver in lines_dict[subf[k-1]].keys():
                        lineloop += [lines_dict[subf[k-1]][ver]+1]
                    elif subf[k-1] in lines_dict[ver].keys():
                        # reversed
                        lineloop += [(lines_dict[ver][subf[k-1]]+1)]
                    else:
                        l = Line([mesh.points[subf[k-1]+1], mesh.points[ver+1]])
                        mesh.addEntity(l)
                        lineloop += [l.nb]
                ll = LineLoop(mesh.getLinesFromIndex(lineloop))
                mesh.addEntity(ll)
                lineloops += [ll.nb]
            s = PlaneSurface([mesh.lineloops[loop] for loop in lineloops])
            mesh.addEntity(s)
            g = mesh.groups.get(self.facetFlags[i])
            if g:
                g.addEntity(s)

    for i, V in enumerate(self.volumes):
        surface_loops = []
        for j, sV in enumerate(V):
            sl = SurfaceLoop((np.array(sV)+1).tolist())
            mesh.addEntity(sl)
            surface_loops += [sl.nb]
        vol = Volume(surface_loops)
        mesh.addEntity(vol)
        g = mesh.groups.get(self.regionFlags[i])
        if g:
            g.addEntity(vol)

    return mesh
