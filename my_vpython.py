from vpython import *
import time


def start_visualization():
    scene.range = 1.2
    scene.caption = 'Right click + drag to rotate camera, scroll wheel to zoom'
    scene.width = 1080
    scene.height = 1080
    scene.ambient = vector(0.4, 0.4, 0.4)

    vs = []
    es = []

    def callback(vertices, edges, vertex_coloring=None):
        if len(vertices) > len(vs):
            for i in range(len(vertices) - len(vs)):
                vs.append(sphere())
        if len(edges) > len(es):
            for i in range(len(edges) - len(es)):
                es.append(cylinder())

        # by default y is locked to be straight up, but I like to use z for that
        for i, v in enumerate(vs):
            if i >= len(vertices):
                v.visible = False
                continue
            vertex = vertices[i]

            r = 0.03
            x, z, y = vertex
            pos = vector(x, y, z)

            v.pos = pos
            v.radius = r
            if vertex_coloring is None:
                v.color = vector(1,1,1)
            else:
                c = vertex_coloring[i]
                v.color = vector(c, 0, 1-c)

        for i, e in enumerate(es):
            if i >= len(edges):
                e.visible = False
                continue
            edge = edges[i]

            ax, az, ay = edge[0]
            bx, bz, by = edge[1]
            dist = sum((edge[1]-edge[0])**2)**0.5

            e.pos = vector(ax, ay, az)
            e.axis = vector(bx-ax, by-ay, bz-az)
            e.length = dist
            e.radius = 0.01
            e.color = vector(0,1,0)

    return callback