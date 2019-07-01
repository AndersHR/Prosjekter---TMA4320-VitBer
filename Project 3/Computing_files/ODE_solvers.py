from numpy import empty, zeros, sqrt, array, float64, copy, vdot, eye, append, arange, swapaxes
from numpy.linalg import solve

#r_0 = [[x_1,y_1],[x_2,y_2], ... ]
#v_0 = [[v_x1,v_y1],[v_x2,v_y2], ... ]
def RK4_exp(n_objects, t_0, t_end, r_0, v_0, m, h, der_r, der_v):
    t = arange(t_0, t_end+h/2, h)
    v_x = empty((n_objects, t.shape[0]))
    v_y = empty((n_objects, t.shape[0]))
    x   = empty((n_objects, t.shape[0]))
    y   = empty((n_objects, t.shape[0]))
    for i in range(n_objects):
        v_x[i, 0] = v_0[i, 0]
        v_y[i, 0] = v_0[i, 1]
        x[i, 0]   = r_0[i, 0]
        y[i, 0]   = r_0[i, 1]
    for i in range(1, t.shape[0]):
        print (i*h)
        x[:, i], y[:, i], v_x[:, i], v_y[:, i] = RK4_exp_step(n_objects, v_x[:, i-1], v_y[:, i-1], x[:, i-1], y[:, i-1], m, h, der_r, der_v)
    return t, x, y, v_x, v_y

def RK4_exp_step(n, v_x, v_y, x, y, m, h, der_r, der_v):
    s1_x, s1_y, s1_vx, s1_vy = empty(n), empty(n), zeros(n), zeros(n)
    s2_x, s2_y, s2_vx, s2_vy = empty(n), empty(n), zeros(n), zeros(n)
    s3_x, s3_y, s3_vx, s3_vy = empty(n), empty(n), zeros(n), zeros(n)
    s4_x, s4_y, s4_vx, s4_vy = empty(n), empty(n), zeros(n), zeros(n)

    #-----------------------------------1-----------------------------------
    for i in range(n):
        s1_x[i] = der_r(v_x[i])
        s1_y[i] = der_r(v_y[i])
        for j in range(i + 1, n):
            a = der_v(x[j] - x[i], y[j] - y[i])
            s1_vx[i] += a * m[j]
            s1_vx[j] -= a * m[i]

            a = der_v(y[j] - y[i], x[j] - x[i])
            s1_vy[i] += a * m[j]
            s1_vy[j] -= a * m[i]

    #-----------------------------------2-----------------------------------
    for i in range(n):
        s2_x[i] = der_r(v_x[i] + (h / 2) * s1_vx[i])
        s2_y[i] = der_r(v_y[i] + (h / 2) * s1_vy[i])
        for j in range(i + 1, n):
            a = der_v(x[j] + (h / 2) * s1_x[j] - (x[i] + (h / 2) * s1_x[i]),
                      y[j] + (h / 2) * s1_y[j] - (y[i] + (h / 2) * s1_y[i]))
            s2_vx[i] += a * m[j]
            s2_vx[j] -= a * m[i]

            a = der_v(y[j] + (h / 2) * s1_y[j] - (y[i] + (h / 2) * s1_y[i]),
                      x[j] + (h / 2) * s1_x[j] - (x[i] + (h / 2) * s1_x[i]))
            s2_vy[i] += a * m[j]
            s2_vy[j] -= a * m[i]

    #-----------------------------------3-----------------------------------
    for i in range(n):
        s3_x[i] = der_r(v_x[i] + (h / 2) * s2_vx[i])
        s3_y[i] = der_r(v_y[i] + (h / 2) * s2_vy[i])

        for j in range(i + 1, n):
            a = der_v(x[j] + (h / 2) * s2_x[j] - (x[i] + (h / 2) * s2_x[i]),
                      y[j] + (h / 2) * s2_y[j] - (y[i] + (h / 2) * s2_y[i]))
            s3_vx[i] += a * m[j]
            s3_vx[j] -= a * m[i]

            a = der_v(y[j] + (h / 2) * s2_y[j] - (y[i] + (h / 2) * s2_y[i]),
                      x[j] + (h / 2) * s2_x[j] - (x[i] + (h / 2) * s2_x[i]))
            s3_vy[i] += a * m[j]
            s3_vy[j] -= a * m[i]

    #-----------------------------------4-----------------------------------
    for i in range(n):
        s4_x[i] = der_r(v_x[i] + h * s3_vx[i])
        s4_y[i] = der_r(v_y[i] + h * s3_vy[i])
        for j in range(i + 1, n):
            a = der_v(x[j] + h * s3_x[j] - (x[i] + h * s3_x[i]),
                      y[j] + h * s3_y[j] - (y[i] + h * s3_y[i]))
            s4_vx[i] += a * m[j]
            s4_vx[j] -= a * m[i]

            a = der_v(y[j] + h * s3_y[j] - (y[i] + h * s3_y[i]),
                      x[j] + h * s3_x[j] - (x[i] + h * s3_x[i]))
            s4_vy[i] += a * m[j]
            s4_vy[j] -= a * m[i]

    return x + (h/6) * (s1_x + 2 * s2_x + 2 * s3_x + s4_x), y + (h/6) * (s1_y + 2 * s2_y + 2 * s3_y + s4_y), \
            v_x + (h / 6) * (s1_vx + 2 * s2_vx + 2 * s3_vx + s4_vx), v_y + (h/6) * (s1_vy + 2 * s2_vy + 2 * s3_vy + s4_vy)

def LG4_imp(n_objects, t_0, t_end, r_0, v_0, m, h, der_r, der_v, der_v_der_primary, der_v_der_secondary):
    MAX_ITER = 10000
    TOL_NEWTON = 1e-10
    c = array([1 / 2 - sqrt(3) / 6, 1 / 2 + sqrt(3) / 6])
    A = array([[1 / 4, 1 / 4 - sqrt(3) / 6], [1 / 4 + sqrt(3) / 6, 1 / 4]])
    b = array([1 / 2, 1 / 2])

    t = array([t_0,], dtype=float64)
    v = empty((1, n_objects, 2), dtype=float64) #lagrer [[vx1, vy1], [vx2, vy2], ....., [vxn, vyn]] per tidsteg
    r = empty((1, n_objects, 2), dtype=float64) #lagrer [a[x1, y1],[x2, y2], .... ,[x_n, y_n]] per tidssteg
    v[0] = v_0
    r[0] = r_0
    while t[-1]<t_end:
        print (t[-1]*1000, "\n")
        #make initial guess
        k_v = zeros(4 * n_objects) #[k1_vx_1, k1_vx_2, k1_vy_1, k1_vy_2, k2_vx_1, k2_vx_2, .....]
        for i in range(n_objects):
            for j in range(i+1, n_objects):
                der_v_x_1 = der_v(r[-1, j, 0] + h*v[-1, j, 0]*c[0] - r[-1, i, 0] - h*v[-1, i, 0]*c[0],
                                  r[-1, j, 1] + h*v[-1, j, 1]*c[0] - r[-1, i, 1] - h*v[-1, i, 1]*c[0])
                der_v_x_2 = der_v(r[-1, j, 0] + h*v[-1, j, 0]*c[1] - r[-1, i, 0] - h*v[-1, i, 0]*c[1],
                                  r[-1, j, 1] + h*v[-1, j, 1]*c[1] - r[-1, i, 1] - h*v[-1, i, 1]*c[1])
                der_v_y_1 = der_v(r[-1, j, 1] + h*v[-1, j, 1]*c[0] - r[-1, i, 1] - h*v[-1, i, 1]*c[0],
                                  r[-1, j, 0] + h*v[-1, j, 0]*c[0] - r[-1, i, 0] - h*v[-1, i, 0]*c[0])
                der_v_y_2 = der_v(r[-1, j, 1] + h*v[-1, j, 1]*c[1] - r[-1, i, 1] - h*v[-1, i, 1]*c[1],
                                  r[-1, j, 0] + h*v[-1, j, 0]*c[1] - r[-1, i, 0] - h*v[-1, i, 0]*c[1])

                k_v[4 * i]     += der_v_x_1 * m[j]
                k_v[4 * i + 1] += der_v_x_2 * m[j]
                k_v[4 * i + 2] += der_v_y_1 * m[j]
                k_v[4 * i + 3] += der_v_y_2 * m[j]

                k_v[4 * j]     -= der_v_x_1 * m[i]
                k_v[4 * j + 1] -= der_v_x_2 * m[i]
                k_v[4 * j + 2] -= der_v_y_1 * m[i]
                k_v[4 * j + 3] -= der_v_y_2 * m[i]
        for iteration in range(MAX_ITER):
            print (iteration)
            #Calculate residue
            resid = copy(k_v)
            k_r = empty(4 * n_objects)
            for i in range(n_objects):
                k_r[4 * i]     = v[-1, i, 0] + h * vdot(A[0, :], k_v[4 * i    :4 * i + 2])
                k_r[4 * i + 1] = v[-1, i, 0] + h * vdot(A[1, :], k_v[4 * i    :4 * i + 2])
                k_r[4 * i + 2] = v[-1, i, 1] + h * vdot(A[0, :], k_v[4 * i + 2:4 * i + 4])
                k_r[4 * i + 3] = v[-1, i, 1] + h * vdot(A[1, :], k_v[4 * i + 2:4 * i + 4])
            for i in range(n_objects):
                for j in range(i+1, n_objects):
                    der_v_x_1 = der_v(r[-1, j, 0] - r[-1, i, 0] + h * (vdot(A[0, :], k_r[4 * j    :4 * j + 2]) - vdot(A[0, :], k_r[4 * i    :4 * i + 2])),
                                      r[-1, j, 1] - r[-1, i, 1] + h * (vdot(A[0, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[0, :], k_r[4 * i + 2:4 * i + 4])))
                    der_v_x_2 = der_v(r[-1, j, 0] - r[-1, i, 0] + h * (vdot(A[1, :], k_r[4 * j    :4 * j + 2]) - vdot(A[1, :], k_r[4 * i    :4 * i + 2])),
                                      r[-1, j, 1] - r[-1, i, 1] + h * (vdot(A[1, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[1, :], k_r[4 * i + 2:4 * i + 4])))
                    der_v_y_1 = der_v(r[-1, j, 1] - r[-1, i, 1] + h * (vdot(A[0, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[0, :], k_r[4 * i + 2:4 * i + 4])),
                                      r[-1, j, 0] - r[-1, i, 0] + h * (vdot(A[0, :], k_r[4 * j    :4 * j + 2]) - vdot(A[0, :], k_r[4 * i    :4 * i + 2])))
                    der_v_y_2 = der_v(r[-1, j, 1] - r[-1, i, 1] + h * (vdot(A[1, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[1, :], k_r[4 * i + 2:4 * i + 4])),
                                      r[-1, j, 0] - r[-1, i, 0] + h * (vdot(A[1, :], k_r[4 * j    :4 * j + 2]) - vdot(A[1, :], k_r[4 * i    :4 * i + 2])))
                    resid[4*i] -= der_v_x_1*m[j]
                    resid[4*i + 1] -= der_v_x_2*m[j]
                    resid[4*i + 2] -= der_v_y_1*m[j]
                    resid[4*i + 3] -= der_v_y_2*m[j]
                    resid[4*j]     += der_v_x_1*m[i]
                    resid[4*j + 1] += der_v_x_2*m[i]
                    resid[4*j + 2] += der_v_y_1*m[i]
                    resid[4*j + 3] += der_v_y_2*m[i]

            #return 1, 1, 1, 1, 1
            #Calculate norm
            max_norm = 0
            for i in range(n_objects):
                x_norm = sqrt(resid[4*i]  **2+resid[4*i+1]**2)/max(abs(v[-1, i, 0]), 1)
                y_norm = sqrt(resid[4*i+2]**2+resid[4*i+3]**2)/max(abs(v[-1, i, 1]), 1)
                if (max_norm < x_norm):
                    max_norm = x_norm
                if (max_norm < y_norm):
                    max_norm = y_norm
            if (max_norm < TOL_NEWTON):
                break

            #Calculate Jacobi
            jac = eye(4 * n_objects)
            for i in range(0, n_objects):
                for j in range(i+1, n_objects):
                    # der f_i1/der_xj
                    val = der_v_der_primary(r[-1, j, 0] - r[-1, i, 0] + h*(vdot(A[0, :], k_r[4 * j    :4 * j + 2]) - vdot(A[0, :], k_r[4 * i    :4 * i + 2])),
                                            r[-1, j, 1] - r[-1, i, 1] + h*(vdot(A[0, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[0, :], k_r[4 * i + 2:4 * i + 4])))
                    jac[4*i, 4*j]   = -h**2*val*(A[0, 0]*A[0, 0] + A[0, 1]*A[1, 0])
                    jac[4*i, 4*j+1] = -h**2*val*(A[0, 0]*A[0, 1] + A[0, 1]*A[1, 1])
                    jac[4*j, 4*i]   = jac[4*i, 4*j]  *m[i]
                    jac[4*j, 4*i+1] = jac[4*i, 4*j+1]*m[i]

                    jac[4*i, 4*i]   -= jac[4*i, 4*j]  *m[j]
                    jac[4*i, 4*i+1] -= jac[4*i, 4*j+1]*m[j]
                    jac[4*j, 4*j]   -= jac[4*i, 4*j]  *m[i]
                    jac[4*j, 4*j+1] -= jac[4*i, 4*j+1]*m[i]

                    jac[4*i, 4*j]   = jac[4*i, 4*j]  *m[j]
                    jac[4*i, 4*j+1] = jac[4*i, 4*j+1]*m[j]

                    # der f_i2/der_xj
                    val = der_v_der_primary(r[-1, j, 0] - r[-1, i, 0] + h*(vdot(A[1, :], k_r[4 * j    :4 * j + 2]) - vdot(A[1, :], k_r[4 * i    :4 * i + 2])),
                                            r[-1, j, 1] - r[-1, i, 1] + h*(vdot(A[1, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[1, :], k_r[4 * i + 2:4 * i + 4])))
                    jac[4*i+1, 4*j]   = -h**2*val*(A[1, 0]*A[0, 0] + A[1, 1]*A[1, 0])
                    jac[4*i+1, 4*j+1] = -h**2*val*(A[1, 0]*A[0, 1] + A[1, 1]*A[1, 1])
                    jac[4*j+1, 4*i]   = jac[4*i, 4*j]  *m[i]
                    jac[4*j+1, 4*i+1] = jac[4*i, 4*j+1]*m[i]

                    jac[4*i+1, 4*i]   -= jac[4*i, 4*j]  *m[j]
                    jac[4*i+1, 4*i+1] -= jac[4*i, 4*j+1]*m[j]
                    jac[4*j+1, 4*j]   -= jac[4*i, 4*j]  *m[i]
                    jac[4*j+1, 4*j+1] -= jac[4*i, 4*j+1]*m[i]

                    jac[4*i+1, 4*j]   = jac[4*i, 4*j]  *m[j]
                    jac[4*i+1, 4*j+1] = jac[4*i, 4*j+1]*m[j]

                    # der f_i3/der_yj
                    val = der_v_der_primary(r[-1, j, 1] - r[-1, i, 1] + h*(vdot(A[0, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[0, :], k_r[4 * i + 2:4 * i + 4])),
                                            r[-1, j, 0] - r[-1, i, 0] + h*(vdot(A[0, :], k_r[4 * j    :4 * j + 2]) - vdot(A[0, :], k_r[4 * i    :4 * i + 2])))
                    jac[4*i+2, 4*j+2] = -h**2*val*(A[0, 0]*A[0, 0] + A[0, 1]*A[1, 0])
                    jac[4*i+2, 4*j+3] = -h**2*val*(A[0, 0]*A[0, 1] + A[0, 1]*A[1, 1])
                    jac[4*j+2, 4*i+2] = jac[4*i+2, 4*j+2]*m[i]
                    jac[4*j+2, 4*i+3] = jac[4*i+2, 4*j+3]*m[i]

                    jac[4*i+2, 4*i+2] -= jac[4*i+2, 4*j+2]*m[j]
                    jac[4*i+2, 4*i+3] -= jac[4*i+2, 4*j+3]*m[j]
                    jac[4*j+2, 4*j+2] -= jac[4*i+2, 4*j+2]*m[i]
                    jac[4*j+2, 4*j+3] -= jac[4*i+2, 4*j+3]*m[i]

                    jac[4*i+2, 4*j+2] = jac[4*i+2, 4*j+2]*m[j]
                    jac[4*i+2, 4*j+3] = jac[4*i+2, 4*j+3]*m[j]

                    # der f_i4/der_yj
                    val = der_v_der_primary(r[-1, j, 1] - r[-1, i, 1] + h*(vdot(A[1, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[1, :], k_r[4 * i + 2:4 * i + 4])),
                                            r[-1, j, 0] - r[-1, i, 0] + h*(vdot(A[1, :], k_r[4 * j    :4 * j + 2]) - vdot(A[1, :], k_r[4 * i    :4 * i + 2])))
                    jac[4*i+3, 4*j+2] = -h**2*val*(A[1, 0]*A[0, 0] + A[1, 1]*A[1, 0])
                    jac[4*i+3, 4*j+3] = -h**2*val*(A[1, 0]*A[0, 1] + A[1, 1]*A[1, 1])
                    jac[4*j+3, 4*i+2] = jac[4*i+3, 4*j+2]*m[i]
                    jac[4*j+3, 4*i+3] = jac[4*i+3, 4*j+3]*m[i]

                    jac[4*i+3, 4*i+2] -= jac[4*i+3, 4*j+2]*m[j]
                    jac[4*i+3, 4*i+3] -= jac[4*i+3, 4*j+3]*m[j]
                    jac[4*j+3, 4*j+2] -= jac[4*i+3, 4*j+2]*m[i]
                    jac[4*j+3, 4*j+3] -= jac[4*i+3, 4*j+3]*m[i]

                    jac[4*i+3, 4*j+2] = jac[4*i+3, 4*j+2]*m[j]
                    jac[4*i+3, 4*j+3] = jac[4*i+3, 4*j+3]*m[j]

                    #----------------------secondary-------------------

                    # der f_i1/der_yj
                    val = der_v_der_secondary(r[-1, j, 0] - r[-1, i, 0] + h*(vdot(A[0, :], k_r[4 * j    :4 * j + 2]) - vdot(A[0, :], k_r[4 * i    :4 * i + 2])),
                                              r[-1, j, 1] - r[-1, i, 1] + h*(vdot(A[0, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[0, :], k_r[4 * i + 2:4 * i + 4])))
                    jac[4*i, 4*j+2] = -h**2*val*(A[0, 0]*A[0, 0] + A[0, 1]*A[1, 0])
                    jac[4*i, 4*j+3] = -h**2*val*(A[0, 0]*A[0, 1] + A[0, 1]*A[1, 1])
                    jac[4*j, 4*i+2] = -jac[4*i, 4*j+2]*m[i]
                    jac[4*j, 4*i+3] = -jac[4*i, 4*j+3]*m[i]

                    jac[4*i, 4*i+2] += jac[4*i, 4*j+2]*m[j]
                    jac[4*i, 4*i+3] += jac[4*i, 4*j+3]*m[j]
                    jac[4*j, 4*j+2] -= jac[4*i, 4*j+2]*m[i]
                    jac[4*j, 4*j+3] -= jac[4*i, 4*j+3]*m[i]

                    jac[4*i, 4*j+2] = jac[4*i, 4*j+2]*m[j]
                    jac[4*i, 4*j+3] = jac[4*i, 4*j+3]*m[j]

                    # der f_i2/der_xj
                    val = der_v_der_secondary(r[-1, j, 0] - r[-1, i, 0] + h*(vdot(A[1, :], k_r[4 * j    :4 * j + 2]) - vdot(A[1, :], k_r[4 * i    :4 * i + 2])),
                                            r[-1, j, 1] - r[-1, i, 1] + h*(vdot(A[1, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[1, :], k_r[4 * i + 2:4 * i + 4])))
                    jac[4*i+1, 4*j+2] = -h**2*val*(A[1, 0]*A[0, 0] + A[1, 1]*A[1, 0])
                    jac[4*i+1, 4*j+3] = -h**2*val*(A[1, 0]*A[0, 1] + A[1, 1]*A[1, 1])
                    jac[4*j+1, 4*i+2] = -jac[4*i, 4*j+2]*m[i]
                    jac[4*j+1, 4*i+3] = -jac[4*i, 4*j+3]*m[i]

                    jac[4*i+1, 4*i+2] += jac[4*i, 4*j+2]*m[j]
                    jac[4*i+1, 4*i+3] += jac[4*i, 4*j+3]*m[j]
                    jac[4*j+1, 4*j+2] -= jac[4*i, 4*j+2]*m[i]
                    jac[4*j+1, 4*j+3] -= jac[4*i, 4*j+3]*m[i]

                    jac[4*i+1, 4*j+2] = jac[4*i, 4*j+2]*m[j]
                    jac[4*i+1, 4*j+3] = jac[4*i, 4*j+3]*m[j]

                    # der f_i3/der_yj
                    val = der_v_der_secondary(r[-1, j, 1] - r[-1, i, 1] + h*(vdot(A[0, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[0, :], k_r[4 * i + 2:4 * i + 4])),
                                            r[-1, j, 0] - r[-1, i, 0] + h*(vdot(A[0, :], k_r[4 * j    :4 * j + 2]) - vdot(A[0, :], k_r[4 * i    :4 * i + 2])))
                    jac[4*i+2, 4*j]   = -h**2*val*(A[0, 0]*A[0, 0] + A[0, 1]*A[1, 0])
                    jac[4*i+2, 4*j+1] = -h**2*val*(A[0, 0]*A[0, 1] + A[0, 1]*A[1, 1])
                    jac[4*j+2, 4*i]   = -jac[4*i+2, 4*j]  *m[i]
                    jac[4*j+2, 4*i+1] = -jac[4*i+2, 4*j+1]*m[i]

                    jac[4*i+2, 4*i]   += jac[4*i+2, 4*j]  *m[j]
                    jac[4*i+2, 4*i+1] += jac[4*i+2, 4*j+1]*m[j]
                    jac[4*j+2, 4*j]   -= jac[4*i+2, 4*j]  *m[i]
                    jac[4*j+2, 4*j+1] -= jac[4*i+2, 4*j+1]*m[i]

                    jac[4*i+2, 4*j]   = jac[4*i+2, 4*j]  *m[j]
                    jac[4*i+2, 4*j+1] = jac[4*i+2, 4*j+1]*m[j]

                    # der f_i4/der_yj
                    val = der_v_der_secondary(r[-1, j, 1] - r[-1, i, 1] + h*(vdot(A[1, :], k_r[4 * j + 2:4 * j + 4]) - vdot(A[1, :], k_r[4 * i + 2:4 * i + 4])),
                                            r[-1, j, 0] - r[-1, i, 0] + h*(vdot(A[1, :], k_r[4 * j    :4 * j + 2]) - vdot(A[1, :], k_r[4 * i    :4 * i + 2])))
                    jac[4*i+3, 4*j]   = -h**2*val*(A[1, 0]*A[0, 0] + A[1, 1]*A[1, 0])
                    jac[4*i+3, 4*j+1] = -h**2*val*(A[1, 0]*A[0, 1] + A[1, 1]*A[1, 1])
                    jac[4*j+3, 4*i]   = -jac[4*i+3, 4*j]  *m[i]
                    jac[4*j+3, 4*i+1] = -jac[4*i+3, 4*j+1]*m[i]

                    jac[4*i+3, 4*i]   += jac[4*i+3, 4*j]  *m[j]
                    jac[4*i+3, 4*i+1] += jac[4*i+3, 4*j+1]*m[j]
                    jac[4*j+3, 4*j]   -= jac[4*i+3, 4*j]  *m[i]
                    jac[4*j+3, 4*j+1] -= jac[4*i+3, 4*j+1]*m[i]

                    jac[4*i+3, 4*j]   = jac[4*i+3, 4*j]  *m[j]
                    jac[4*i+3, 4*j+1] = jac[4*i+3, 4*j+1]*m[j]
            k_v = k_v - solve(jac, resid)

        if iteration==(MAX_ITER-1):
            print ("\n\nHOLY_COW\n\n")

        t = append(t, t[-1]+h)
        v = append(v, [empty((n_objects, 2))], axis=0)
        r = append(r, [empty((n_objects, 2))], axis=0)

        for i in range(n_objects):
            v[-1, i, 0] = v[-2, i, 0] + h*vdot(b, k_v[4*i  : 4*i+2])
            v[-1, i, 1] = v[-2, i, 1] + h*vdot(b, k_v[4*i+2: 4*i+4])
            r[-1, i, 0] = r[-2, i, 0] + h * vdot(b, k_r[4 * i: 4 * i + 2])
            r[-1, i, 1] = r[-2, i, 1] + h * vdot(b, k_r[4 * i + 2: 4 * i + 4])




    return t, swapaxes(r[:, :, 0], 0, 1), swapaxes(r[:, :, 1], 0, 1), swapaxes(v[:, :, 0], 0, 1), swapaxes(v[:, :, 1], 0, 1)

def EC1_exp(n_objects, t_0, t_end, r_0, v_0, m, h, der_r, der_v):
    t = arange(t_0, t_end+h/2, h)
    v = empty((t.shape[0], n_objects, 2), dtype=float64)  # lagrer [[vx1, vy1], [vx2, vy2], ....., [vxn, vyn]] per tidsteg
    r = empty((t.shape[0], n_objects, 2), dtype=float64)  # lagrer [a[x1, y1],[x2, y2], .... ,[x_n, y_n]] per tidssteg
    v[0] = v_0
    r[0] = r_0
    for iter_index in range(1, t.shape[0]):
        print (iter_index)
        acc = zeros((n_objects, 2))
        for i in range(n_objects):
            for j in range(i+1, n_objects):
                a = der_v(r[iter_index-1, j, 0] - r[iter_index-1, i, 0], r[iter_index-1, j, 1] - r[iter_index-1, i, 1])
                acc[i, 0] += a * m[j]
                acc[j, 0] -= a * m[i]
                a = der_v(r[iter_index-1, j, 1] - r[iter_index-1, i, 1], r[iter_index-1, j, 0] - r[iter_index-1, i, 0])
                acc[i, 1] += a * m[j]
                acc[j, 1] -= a * m[i]

        v[iter_index] = array(v[iter_index-1]+h*acc)
        r[iter_index] = array(r[iter_index-1]+h*v[iter_index])

    return t, swapaxes(r[:, :, 0], 0, 1), swapaxes(r[:, :, 1], 0, 1), swapaxes(v[:, :, 0], 0, 1), swapaxes(v[:, :, 1], 0, 1)