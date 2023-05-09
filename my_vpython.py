from vpython import *
import time


def start_visualization(vertices, edges, vertex_coloring=None):
    scene.range = 1.2
    scene.caption = 'Right click + drag to rotate camera, scroll wheel to zoom'
    scene.width = 1080
    scene.height = 1080
    scene.ambient = vector(0.4, 0.4, 0.4)

    # by default y is locked to be straight up, but I like to use z for that
    vs = []
    for i, vertex in enumerate(vertices):
        r = 0.03
        x, z, y = vertex
        pos = vector(x, y, z)
        if vertex_coloring is None:
            v = sphere(pos=pos, radius=r)
        else:
            c = vertex_coloring[i]
            v = sphere(pos=pos, radius=r, color=vector(c, 0, 1-c))
        vs.append(v)

    es = []
    for edge in edges:
        ax, az, ay = edge[0]
        bx, bz, by = edge[1]
        dist = sum((edge[1]-edge[0])**2)**0.5
        e = cylinder(pos=vector(ax, ay, az), axis=vector(bx-ax, by-ay, bz-az), length=dist, radius=0.01, color=vector(0,1,0))
        es.append(e)


    #while True:
    #    #rate(60)
    #    time.sleep(1)

