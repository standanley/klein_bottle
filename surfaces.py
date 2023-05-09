import numpy as np

def surface_sphere(uv):
    u, v = uv
    x = np.cos(2 * np.pi * u) * np.sin(np.pi * v)
    y = np.sin(2 * np.pi * u) * np.sin(np.pi * v)
    z = np.cos(np.pi * v)
    return np.array([x, y, z])

def surface_cylinder(uv):
    r1, r2 = 0.5, 0.2
    u, v = uv
    x = (r1 + r2*np.cos(2 * np.pi * v)) * np.cos(2 * np.pi * u)
    y = (r1 + r2*np.cos(2 * np.pi * v)) * np.sin(2 * np.pi * u)
    z = r2 * np.sin(2 * np.pi * v)
    return np.array([x, y, z])

def surface_klein(uv):
    h_top = 0.5
    h_bottom = 1
    w_right = 0.5
    w_left = 0.3
    r1 = 0.1
    r2 = 0.4
    r3 = 0.6

    # difficulty with this parameterization:
    # consider what happens as the inner tube approaches the bottom:
    # r must be growing as u_z is at a standstill
    # u_z must be growing as u_x is at a standstill
    # so I don't think we can accomplish this with trig functions along
    # solution: we will use nested trig functions for u_x to make it doubly slow

    def jug_handle(t, h_top, h_bottom, w):
        # as t ranges from 0 to 1, this ranges from (0, h_top) to (0, -h_bottom) with a stop at (w, 0)
        # I think h_bottom should equal 2*h_top to have continuous second derivative, and w=h_top looks good
        if t < 0.5:
            t0 = 2 * t
            x = w * np.sin(1 / 2 * np.pi * t0)
            z = h_top * np.cos(1/2 * np.pi * t0)
        else:
            t1 = 2 * t - 1
            # as tu goes from 0 to 1 so does tu_stretch, but it slows at the end
            t_stretch = np.sin(np.pi / 2 * t1)
            x = w * (0.5 * np.cos(np.pi * t_stretch) + 0.5)
            z = h_bottom * np.cos(1/2 * np.pi * (t1+1))
        return x, z

    def partial_f(f, t):
        # numerical derivative of f with respect to t at t
        eps = 1e-6
        t_low = max(0, t-eps)
        t_high = min(1, t+eps)
        x_low, z_low = f(t_low)
        x_high, z_high = f(t_high)
        dx = (x_high - x_low) / (t_high - t_low)
        dz = (z_high - z_low) / (t_high - t_low)
        return dx, dz

    def cylinderify(f, t, v, r):
        x_u, z_u = f(t)
        y_u = 0

        dx, dz = partial_f(f, t)
        ds = np.sqrt(dx**2 + dz**2)

        x_v = r * np.cos(2*np.pi*v)
        y_v = r * np.sin(2*np.pi*v)

        x = x_u + x_v * (dz/ds)
        y = y_u + y_v
        z = z_u + x_v * (-dx/ds)
        return x, y, z


    u, v = uv
    if u < 0.5:
        # right hand jug handle
        t = 2*u
        f = lambda t0: jug_handle(t0, h_top, h_bottom, w_right)

        r = r1 if t < 0.5 else r1 + (1-np.cos(np.pi*(t-0.5)))*(r2 - r1)
        x, y, z = cylinderify(f, t, v, r)

    else:
        # left hand jug handle
        t = 2*u - 1
        f = lambda t0: jug_handle(1-t0, h_top, h_bottom, -w_left)

        r = r2 + np.sin(2*np.pi*t)*(r3-r2) if t < 0.5 else r2 + np.sin(np.pi*(t-0.5))*(r1-r2)
        # because the direction of f reverses at u=0.5, we need to do something to v so that
        # the seams line up correctly here. Note that dx/dt continues to be right-to-left, while
        # dz/dt switches from down to up, so to make the circles consistent we need to change
        # our process for v. I don't know why it is the way it is exactly
        x, y, z = cylinderify(f, t, (-v+0.5), r)


    return np.array([x, y, z])
