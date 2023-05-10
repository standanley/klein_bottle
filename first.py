import itertools
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
    return uvs[:, 1]

def undo_uv_wrapping(uvs):
    # given a set of uvs that are reasonably close together, make sure they aren't
    # across a seam
    # i.e. change 0.1 and 0.9 into 0.1 and -0.1
    # this probably breaks if there are uvs in the -0.1 to 0.1 range and more in the 0.4 to 0.6 range;
    # we assume they have been pre-filtered so they are already reasonably close together
    uvs = np.array(uvs)
    uvs_final = np.copy(uvs)
    u_rep = uvs_final[0][0]
    if u_rep < 0.25 or u_rep > 0.75:
        # kind of near border, let's use border version, everyone must be in [-0.5, 0.5]
        for i in range(len(uvs_final)):
            if uvs_final[i][0] > 0.5:
                uvs_final[i][0] -= 1
                # subtracting 1 always crosses the border, so we need to flip v
                uvs_final[i][1] = (0.5 - uvs_final[i][1])%1

    # same game, except shifting v never forces us to flip u
    # we want to be more careful with the v, and put the split at the widest break
    temp = uvs_final.copy()
    v_rep = uvs_final[0][1]

    vs = sorted(uvs_final[:,1])
    vs_gaps = [(vs[(i+1)%len(vs)]-vs[i])%1 for i in range(len(vs))]
    max_gap_location = np.argmax(vs_gaps)
    break_location = (vs[max_gap_location] + vs_gaps[max_gap_location] / 2) % 1

    for i in range(len(uvs_final)):
        if uvs_final[i][1] > break_location:
            uvs_final[i][1] -= 1
    us = uvs_final[:,0]
    vs = uvs_final[:,1]
    # these assertions can absolutely fail because sometimes the v is that bad
    # in other words, we sometimes get an edge spanning right through the middle of the
    # thin pipe, and then it's hard to say which side it should be counted as going around
    assert 0 < max(us) - min(us) < 0.2
    #assert 0 < max(vs) - min(vs) < 0.4
    # (not (0 < max(vs) - min(vs) < 0.1)) and 0.9>us[0]>0.6
    return uvs_final



def remove_crosses(uvs, xyzs, pairs):
    edge_lookup = {(i, j): dist for i, j, dist in pairs}
    edges_final = set(edge_lookup.keys())
    print('starting with', len(edges_final), 'edges')
    #edges_final = []

    def check_cross(a_i, a_j, b_i, b_j):
        # first we need to undo any uv wrapping
        # unwrap u first, becasue
        uvs_quad = [uvs[k] for k in [a_i, a_j, b_i, b_j]]
        uvs_unwrapped = undo_uv_wrapping(uvs_quad)
        #print(uvs_unwrapped)

        def is_left(point, v, test_point):
            # return true iff test_point is to the left looking from point in direction v
            w = test_point - point
            cross_product = v[0]*w[1] - v[1]*w[0]
            return cross_product < 0

        e, f, g, h = uvs_unwrapped
        vec_a = f - e
        vec_b = h - g
        cross = ((is_left(e, vec_a, g) != is_left(e, vec_a, h))
                 and (is_left(g, vec_b, e) != is_left(g, vec_b, f)))
        return cross


    for a in range(len(pairs) - 1):
        for b in range(a+1, len(pairs)):
            # we rely on the fact that if a and b cross, the surrounding square
            # all has shorter edges than the crossing diagonals
            # so we can quickly rule out pairs of edges if the surrounding square isnt in the lookup
            a_i, a_j, a_dist = pairs[a]
            b_i, b_j, b_dist = pairs[b]

            if (a_i, a_j) not in edges_final or (b_i, b_j) not in edges_final:
                # one of these has already been removed for crossing something else
                continue

            can_skip = False
            for test1, test2 in itertools.product((a_i, a_j), (b_i, b_j)):
                if tuple(sorted((test1, test2))) not in edge_lookup:
                    can_skip = True
                    break
            if can_skip:
                continue

            if not check_cross(a_i, a_j, b_i, b_j):
                continue

            if a_dist < b_dist:
                edges_final.remove((b_i, b_j))
            else:
                edges_final.remove((a_i, a_j))
            #edges_final.append((a, b))

    print('in the end, ', len(edges_final), 'edges left')
    return list(edges_final)

N = 200
uvs = np.random.uniform(size=(N, 2))

spring_max_length = 0.3
spring_k = 0.05


callback = start_visualization()
xyzs = surface(uvs)
callback(xyzs, [], get_colors(uvs))

#for i in range(200):
#    jac = jacobian(uvs)
#    xyzs = surface(uvs)
#    print('starting pairs')
#    pairs = get_pairs(uvs, xyzs, spring_max_length)
#    print('finished pairs')
#    edges = [(xyzs[i], xyzs[j]) for i, j, dist in pairs]
#    callback(xyzs, edges, get_colors(uvs))
#    if i == 0:
#        print('sleeping 10 seconds to let vpython load better')
#        time.sleep(10)
#
#    uvs = spring_push(uvs, xyzs, jac, pairs, spring_max_length, spring_k)

#np.savetxt('uvs_moved.csv', uvs, delimiter=',')
uvs = np.loadtxt('uvs_moved.csv', delimiter=',')

xyzs = surface(uvs)
pairs = get_pairs(uvs, xyzs, spring_max_length)
edges = [(xyzs[i], xyzs[j]) for i, j, dist in pairs]
callback(xyzs, edges, get_colors(uvs))

edges_final = remove_crosses(uvs, xyzs, pairs)
colors_temp = np.zeros(len(pairs))
for a, b in edges_final:
    colors_temp[a] = 1
    colors_temp[b] = 1

vertex_coloring = get_colors(uvs)
#callback(xyzs, edges, vertex_coloring, colors_temp)

edges2 = [(xyzs[a], xyzs[b]) for a, b in edges_final]
callback(xyzs, edges2, vertex_coloring)
