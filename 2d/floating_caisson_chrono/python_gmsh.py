
filename = 'test1.msh'

mshfile = open(filename, 'r')
nodes = []
edges = []
triangles = []
tetrahedra = []
triangle_nb = 0
edge_nb = 0

nd = 2

for line in mshfile:
    switch = None
    if 'Nodes' in line:
        switch = 'nodes'
        switch_count = 0
    if 'Elements' in line:
        switch = 'elements'
        switch_count = 0
    if switch == 'nodes':
        words = line.split()
        if switch_count == 0:
            node_nb = float(words[0])
        else:
            nid = int(words[0])
            x, y, z = float(words[1]), float(words[2]), float(words[3])
            nodes += [[node_nb, x, y, z]]
    if switch == 'elements':
        words = words.split()
        if switch_count == 0:
            el_nb = float(words[0])
        else:
            el_id = int(words[0])
            el_type = int(words[1])
            nb_tags = int(words[2])
            if nb_tags == 2:
                flag = int(words[3])
            else:
                flag = 0
            s = 3+nb_tags # starting index on words for element info
            if el_type == 1: # segment
                edge_nb += 1
                edges += [[edge_nb, int(words[s]), int(words[s+1]), flag]]
            elif el_type == 2: # triangle
                triangle_nb += 1
                triangles += [[triangle_nb, int(words[s]), int(words[s+1]), flag]]
            # elif el_type == 15: # node
mshfile.close()
nodes = np.array(nodes)


# with open(filename, 'w') as elefile
#     writer = csv.writer()

header = '{0:d} {1:d} 0 1'.format(node_nb, nd)
np.savetxt('mesh.edge', nodes[:,:nd+1], fmt='%d', header=header, comments='')

header = '{0:d} 1'.format(edge_nb)
np.savetxt('mesh.edge', edges, fmt='%d', header=header, comments='')

if nd == 2:
    header = '{0:d} 3 1'.format(triangle_nb)
    np.savetxt('mesh.ele', triangles, fmt='%d', header=header, comments='')
if nd == 3:
    header = '{0:d} 4 1'.format(tretra_nb)
    np.savetxt('mesh.ele', tetrahedra, fmt='%d', header=header, comments='')
    
            
            



    switch_count += 1
