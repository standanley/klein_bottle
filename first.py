import time

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
    # return xyz distance between things; or float(inf) if u distance is greater than u_limit
    # beware: if you want to add a v_limit, you need to think about the seams, and it's difficult
    u_limit = 0.15
    u_a, u_b = uv_a[0], uv_b[0]
    if u_limit < abs(u_b - u_a) < 1-u_limit:
        return INF

    return np.sqrt(np.sum((xyz_b - xyz_a)**2))

def get_pairs(uvs, xyzs, maximum):
    close = []
    for i in range(len(uvs)-1):
        for j in range(i+1, len(uvs)):
            uv_i, uv_j = uvs[i], uvs[j]
            xyz_i, xyz_j = xyzs[i], xyzs[j]
            dist = get_dist(uv_i, uv_j, xyz_i, xyz_j)
            if dist < maximum:
                close.append((i, j, dist))
    #print('found', len(close), 'edges')
    return close

    #for i, j, dist in close:
    #    xyz_i, xyz_j = xyzs[i], xyzs[j]

def bump_uv(uv, duv):
    # this would be simple but for wrapping behavior, particularly because of the non-orientatableness
    u_temp = (uv[0] + duv[0]) % 2
    u_final = u_temp % 1

    v_temp = (uv[1] + duv[1]) % 1
    v_final = (0.5 - v_temp)%1 if u_temp >= 1 else v_temp

    return np.array([u_final, v_final])


def jacobian(uvs):
    # input is Nx2
    # output is Nx3x2
    eps = 1e-6

    def _jacobian(uv):
        ans = np.zeros((3,2))
        for i in range(2):
            bump_low = [0, 0]
            bump_low[i] = -eps
            bump_high = [0, 0]
            bump_high[i] = eps
            xyz_low = surface(bump_uv(uv, bump_low))
            xyz_high = surface(bump_uv(uv, bump_high))
            ans[:,i] = (xyz_high - xyz_low)/(2*eps)
        return ans

    ans = np.vectorize(_jacobian, signature='(2)->(3,2)')(uvs)
    return ans


def spring_push(uvs, xyzs, jac, edges, max_dist, k):
    # we don't do any momentum, so there's no timestep or mass, just dx=f
    forces_felt = np.zeros((len(uvs), 3))
    spring_count = 0
    total_dist = 0
    for i, j, dist in edges:
        if dist > max_dist:
            continue
        spring_count += 1
        total_dist += dist
        f_magnitude = k*(max_dist - dist)
        # points from i to j
        f_vec = (xyzs[j] - xyzs[i])/dist * f_magnitude
        forces_felt[i] += -f_vec
        forces_felt[j] += f_vec

    print('spring count:', spring_count, '\tavg spring length', total_dist/spring_count)

    uvs_new = np.zeros(uvs.shape)
    for i in range(len(uvs)):
        if np.sum(forces_felt**2) == 0:
            uvs_new[i, :] = uvs[i, :]
        else:
            # no momentum, just move it by the force vector
            f_xyz = forces_felt[i]
            assert f_xyz.shape == (3,)
            inverse_jac= np.linalg.pinv(jac[i, :, :])
            assert inverse_jac.shape == (2, 3)
            f_uv = inverse_jac @ f_xyz
            assert f_uv.shape == (2,)
            uv_new = bump_uv(uvs[i], f_uv)
            uvs_new[i, :] = uv_new

    return uvs_new

def get_colors(uvs):
    # use channel 1 to see the seam
    return uvs[:, 0]

N = 200
uvs = np.random.uniform(size=(N, 2))

spring_max_length = 0.3
spring_k = 0.05


callback = start_visualization()
xyzs = surface(uvs)
callback(xyzs, [], get_colors(uvs))

for i in range(200):
    jac = jacobian(uvs)
    xyzs = surface(uvs)
    print('starting pairs')
    pairs = get_pairs(uvs, xyzs, spring_max_length)
    print('finished pairs')
    edges = [(xyzs[i], xyzs[j]) for i, j, dist in pairs]
    callback(xyzs, edges, get_colors(uvs))
    if i == 0:
        print('sleeping 10 seconds to let vpython load better')
        time.sleep(10)

    uvs = spring_push(uvs, xyzs, jac, pairs, spring_max_length, spring_k)

xyzs = surface(uvs)
pairs = get_pairs(uvs, xyzs, spring_max_length)
edges = [(xyzs[i], xyzs[j]) for i, j, dist in pairs]
callback(xyzs, edges, get_colors(uvs))

np.savetxt('uvs_moved.csv', uvs, delimiter=',')
uvs = np.loadtxt('uvs_moved.csv', delimiter=',')

