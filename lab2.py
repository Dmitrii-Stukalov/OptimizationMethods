import numpy as np


def func1(x):
    return (-5 * x ** 5) + (4 * x ** 4) - (12 * x ** 3) + (11 * x ** 2) - (2 * x) + 1


def func1_der(x):
    return -25 * x ** 4 + 16 * x ** 3 - 36 * x ** 2 + 22 * x - 2


def func2(x):
    return -3 * x * np.sin(0.75 * x) + np.exp(-2 * x)


def func2_der(x):
    return -3 * (np.sin(0.75 * x) + 0.75 * np.cos(0.75 * x) * x) - 2 * np.exp(-2 * x)


def func3(x):
    return np.exp(3 * x) + 5 * np.exp(-2 * x)


def func3_der(x):
    return 3 * np.exp(3 * x) - 10 * np.exp(-2 * x)


def secant(x, w, f_x, f_w):
    return (x * f_w - w * f_x) / (f_w - f_x)


def brent(a, c, eps, f, f_der):
    x = w = v = (a + c) / 2
    f_x = f_w = f_v = f(x)
    f_der_x = f_der_w = f_der_v = f_der(x)
    d = e = c - a
    u = None
    n = 0
    while abs(c - a) > eps:
        n += 1
        g = e
        e = d
        prev_u = u
        u = None
        if x != w and f_der_x != f_der_w:
            u1 = secant(x, w, f_x, f_w)

            if a + eps <= u1 <= c - eps and abs(u1 - x) < g / 2:
                u = u1
        if x != v and f_der_x != f_der_v:
            u2 = secant(x, v, f_x, f_v)

            if a + eps <= u2 <= c - eps and abs(u2 - x) < g / 2:
                if u is None or abs(u2 - x) < abs(u - x):
                    u = u2
        if u is None:
            if f_der_x > 0:
                u = (a + x) / 2
            else:
                u = (x + c) / 2
        if abs(u - x) < eps:
            u = x + np.sign(u - x) * eps

        d = abs(x - u)
        f_u = f(u)
        f_der_u = f_der(u)
        if f_u <= f_x:
            if u >= x:
                a = x
            else:
                c = x
            v = w
            w = x
            x = u
            f_v = f_w
            f_w = f_x
            f_x = f_u
            f_der_v = f_der_w
            f_der_w = f_der_x
            f_der_x = f_der_u
        else:
            if u >= x:
                c = u
            else:
                a = u
            if f_u <= f_w or w == x:
                v = w
                w = u
                f_v = f_w
                f_w = f_u
                f_der_v = f_der_w
                f_der_w = f_der_u
            elif f_u <= f_v or v == x or v == w:
                v = u
                f_v = f_u
                f_der_v = f_der_u
        if n > 1 and abs(prev_u - u) < eps:
            break
    return (a + c) / 2

# a, b = -0.5, 0.5
# e = 1e-5
# print(brent(a, b, e, func1, func1_der))

# a, b = 0, 2 * np.pi
# e = 1e-5
# print(brent(a, b, e, func2, func2_der))

# a, b = 0, 1
# e = 1e-5
# print(brent(a, b, e, func3, func3_der))


def f1(x1, x2):
    return 100 * (x2 - x1 ** 2) ** 2 + (1 - x1) ** 2


def f1_der(x1, x2, index):
    return 2 * (200 * x1 ** 3 - 200 * x1 * x2 + x1 - 1) if index == 1 else 200 * (x2 - x1 ** 2)


def f2(x1, x2):
    return (x2 - x1 ** 2) ** 2 + (1 - x1) ** 2


def f2_der(x1, x2, index):
    return 2 * (2 * x1 ** 3 - 2 * x1 * x2 + x1 - 1) if index == 1 else 2 * (x2 - x1 ** 2)


def f3(x1, x2):
    return (1.5 - x1 * (1 - x2)) ** 2 + (2.25 - x1 * (1 - x2 ** 2)) ** 2 + (2.625 - x1 * (1 - x2 ** 3)) ** 2


def f3_der(x1, x2, index):
    return -12.75 + 3 * x2 + 4.5 * x2 ** 2 + 5.25 * x2 ** 3 + 2 * x1 * (
            3 - 2 * x2 - x2 ** 2 - 2 * x2 ** 3 + x2 ** 4 + x2 ** 6) if index == 1 else x1 * (
            3 + 9 * x2 + 15.75 * x2 ** 2 + x1 * (-2 - 2 * x2 - 6 * x2 ** 2 + 4 * x2 ** 3 + 6 * x2 ** 5))


def f4(x1, x2, x3, x4):
    return (x1 + x2) ** 2 + 5 * (x3 - x4) ** 2 + (x2 - 2 * x3) ** 4 + 10 * (x1 - x4) ** 4


def f4_der(x1, x2, x3, x4, index):
    if index == 1:
        return 2 * (20 * (x1 - x4) ** 3 + x1 + x2)
    elif index == 2:
        return 2 * (2 * (x2 - 2 * x3) ** 3 + x1 + x2)
    elif index == 3:
        return 10 * (x3 - x4) - 8 * (x2 - 2 * x3) ** 3
    else:
        return 10 * (-x3 - 4 * (x1 - x4) ** 3 + x4)


def coordinate_descent(f):
    x = []
    return x

# fig = plt.figure()
# ax = fig.gca(projection='3d')

# X = np.arange(0, 2, 0.01)
# Y = np.arange(0, 2, 0.01)
# X, Y = np.meshgrid(X, Y)
# Z = f1(X, Y)

# surf = ax.plot_surface(X, Y, Z, cmap=cm.get_cmap('coolwarm'),
#                        linewidth=0, antialiased=False)
# ax.view_init(60, 35)
# plt.show()

