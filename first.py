import numpy as np
np.random.seed(4)
from my_vpython import start_visualization
from surfaces import surface_klein

INF = float('inf')

print('hi')

#start_visualization([
#    [0, 0, 0],
#    [0.5, 0, 0],
#    [0, 1, 0]
#])

#N = 50
#us, vs = np.meshgrid(*[np.arange(0, 1, 1/N)]*2)
#uvs = np.stack((us, vs), 2).reshape((-1, 2))
surface = np.vectorize(surface_klein, signature='(2)->(3)')

def get_dist(uv_a, uv_b, xyz_a, xyz_b):
    # return xyz distance between things; or float(inf) if u distance is greater than 0.25
    u_a, u_b = uv_a[0], uv_b[0]
    if 0.25 < abs(u_b - u_a) < 0.75:
        return INF

    return np.sqrt(np.sum((xyz_b - xyz_a)**2))

def get_pairs(uvs):
    xyzs = surface(uvs)
    close = []
    maximum = 0.15
    for i in range(len(uvs)-1):
        for j in range(i+1, len(uvs)):
            uv_i, uv_j = uvs[i], uvs[j]
            xyz_i, xyz_j = xyzs[i], xyzs[j]
            dist = get_dist(uv_i, uv_j, xyz_i, xyz_j)
            if dist < maximum:
                close.append((i, j, dist))
    print('found', len(close), 'edges')
    return close

    #for i, j, dist in close:
    #    xyz_i, xyz_j = xyzs[i], xyzs[j]

N = 500
uvs = np.random.uniform(size=(N, 2))

pairs = get_pairs(uvs)
xyzs = surface(uvs)
edges = [(xyzs[i], xyzs[j]) for i, j, dist in pairs]

start_visualization(xyzs, edges)
