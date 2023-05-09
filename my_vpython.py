from vpython import *
import time


def start_visualization(vertices, edges):
    scene.range = 1.2
    scene.caption = 'Right click + drag to rotate camera, scroll wheel to zoom'
    scene.width = 1080
    scene.height = 1080
    scene.ambient = vector(0.4, 0.4, 0.4)

    # by default y is locked to be straight up, but I like to use z for that
    vs = [sphere(pos=vector(x, z, y), radius=0.03) for x, y, z in vertices]
    es = []
    for edge in edges:
        ax, az, ay = edge[0]
        bx, bz, by = edge[1]
        dist = sum((edge[1]-edge[0])**2)**0.5
        e = cylinder(pos=vector(ax, ay, az), axis=vector(bx-ax, by-ay, bz-az), length=dist, radius=0.01, color=vector(0,0,1))
        es.append(e)


    #while True:
    #    #rate(60)
    #    time.sleep(1)

