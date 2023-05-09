from vpython import *
import time


def start_visualization(points):
    scene.range = 1.2
    scene.caption = 'Right click + drag to rotate camera, scroll wheel to zoom'
    scene.width = 1080
    scene.height = 1080
    scene.ambient = vector(0.4, 0.4, 0.4)

    # by default y is locked to be straight up, but I like to use z for that
    spheres = [sphere(pos=vector(x, z, y), radius=0.03) for x, y, z in points]

    #while True:
    #    #rate(60)
    #    time.sleep(1)

