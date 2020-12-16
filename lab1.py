from math import log, sin, exp, sqrt, pi, cos

import matplotlib.pyplot as plt
import numpy as np


def f1(x):
    return (-5 * x ** 5) + (4 * x ** 4) - (12 * x ** 3) + (11 * x ** 2) - (2 * x) + 1


def f2(x):
    return log(x - 2) ** 2 + log(10 - x) ** 2 - x ** 0.2


def f3(x):
    return -3 * x * sin(0.75 * x) + exp(-2 * x)


def f4(x):
    return exp(3 * x) + 5 * exp(-2 * x)


def f5(x):
    return 0.2 * x * log(x) + (x - 2.3) ** 2


def f6(x):
    return -2 * x ** 4 - 3 * x ** 3 + 10 * x ** 2 - x + 5


def dichotomy(a, b, e, f, it=0, prev=1):
    x1 = (a + b) / 2 - e / 3
    x2 = (a + b) / 2 + e / 3
    y1 = f(x1)
    y2 = f(x2)
    # print(it, a, b, prev / (b - a), x1, x2, y1, y2)
    if b - a < e:
        # print(it)
        return (x1 + x2) / 2
    elif y1 > y2:
        return dichotomy(x1, b, e, f, it + 1, b - a)
    elif y1 < y2:
        return dichotomy(a, x2, e, f, it + 1, b - a)
    elif y1 == y2:
        return dichotomy(x1, x2, e, f, it + 1, b - a)


def golden_section(a, b, x1, x2, y1, y2, e, f, iter=0, last=1):
    # print(iter, a, b, (b - a) / last, x1, x2, y1, y2)
    if b - a < e:
        # print(iter)
        return (x1 + x2) / 2
    elif y1 > y2:
        p_a = a
        a = x1
        x1 = x2
        y1 = y2
        x2 = a + 0.618 * (b - a)
        y2 = f(x2)
        return golden_section(a, b, x1, x2, y1, y2, e, f, iter + 1, b - p_a)
    elif y1 <= y2:
        p_b = b
        b = x2
        x2 = x1
        y2 = y1
        x1 = a + 0.381 * (b - a)
        y1 = f(x1)
        return golden_section(a, b, x1, x2, y1, y2, e, f, iter + 1, p_b - a)


def fib(n):
    n += 1
    return (((1 + sqrt(5)) / 2) ** n - ((1 - sqrt(5)) / 2) ** n) / sqrt(5)


def fibonacci(a, b, e, f):
    n = 0
    while fib(n) <= (b - a) / e:
        n += 1

    x1 = a + fib(n - 2) / fib(n) * (b - a)
    x2 = a + fib(n - 1) / fib(n) * (b - a)
    y1 = f(x1)
    y2 = f(x2)
    prev = 1
    for k in range(1, n - 1):
        # print(k - 1, a, b, (b - a) / prev, x1, x2, y1, y2)
        prev = b - a
        if y1 > y2:
            a = x1
            x1 = x2
            y1 = y2
            x2 = a + fib(n - k - 1) / fib(n - k) * (b - a)
            y2 = f(x2)
        elif y1 <= y2:
            b = x2
            x2 = x1
            y2 = y1
            x1 = a + fib(n - k - 2) / fib(n - k) * (b - a)
            y1 = f(x1)
    return (x1 + x2) / 2
    # return n - 3


def parabola(x1, x2, x3, e, f, it=0, prev=1):
    y1 = f(x1)
    y2 = f(x2)
    y3 = f(x3)

    u = parabola_min(x1, x2, x3, y1, y2, y3)
    y_u = f(u)
    # print(it, x1, x3, (x3 - x1) / prev, x2, u, y2, y_u)

    if abs(x2 - u) < e:
        # print(it)
        return (u + x2) / 2
    elif y2 >= y_u:
        if u >= x2:
            return parabola(x2, u, b, e, f, it + 1, x3 - x1)
        else:
            return parabola(a, u, x2, e, f, it + 1, x3 - x1)
    elif y2 < y_u:
        if u >= x2:
            return parabola(a, x2, u, e, f, it + 1, x3 - x1)
        else:
            return parabola(u, x2, b, e, f, it + 1, x3 - x1)


def parabola_min(x1, x2, x3, y1, y2, y3):
    try:
        return x2 - ((x2 - x1) ** 2 * (y2 - y3) - (x2 - x3) ** 2 * (y2 - y1)) / (
                2 * ((x2 - x1) * (y2 - y3) - (x2 - x3) * (y2 - y1)))
    except:
        return x2


def brent(a, c, eps, f):
    K = (3 - sqrt(5)) / 2
    x = w = v = (a + c) / 2
    f_x = f_w = f_v = f(x)
    d = e = c - a
    u = None
    n = 0
    prev = 1
    while abs(c - a) > eps:
        # print(n, a, c, (c - a) / prev, x, w, v, u, f_x, f_w, f_v, f(u) if u is not None else '-')
        n += 1
        g = e
        e = d
        prev_u = u
        if x != w and w != v and x != v and f_x != f_w and f_x != f_v and f_w != f_v:
            if x < w < v:
                u = parabola_min(x, w, v, f_x, f_w, f_v)
            elif v < w < x:
                u = parabola_min(v, w, x, f_v, f_w, f_x)
            elif x < v < w:
                u = parabola_min(x, v, w, f_x, f_v, f_w)
            elif v < x < w:
                u = parabola_min(v, x, w, f_v, f_x, f_w)
            elif w < v < x:
                u = parabola_min(w, v, x, f_w, f_v, f_x)
            elif w < x < v:
                u = parabola_min(w, x, v, f_w, f_x, f_v)

            if a + eps <= u <= c - eps and abs(u - x) < g / 2:
                d = abs(u - x)
            else:
                u = None
        else:
            u = None
        if u is None:
            if x < (c + a) / 2:
                u = x + K * (c - x)
                d = c - x
            else:
                u = x - K * (x - a)
                d = x - a
        if abs(u - x) < eps:
            u = x + np.sign(u - x) * eps

        f_u = f(u)
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
            elif f_u <= f_v or v == x or v == w:
                v = u
                f_v = f_u
        if n > 1 and abs(prev_u - u) < eps:
            break
    # print(n)
    return (a + c) / 2


def draw_plot(x, y):
    plt.plot(x, y)
    plt.show()


def test_func(a, b, e, f):
    x = np.arange(a, b, e)
    y = [f(x_i) for x_i in x]
    draw_plot(x, y)

    print(f.__name__)

    print('Ideal:', x[list(y).index(min(y))])

    print('Dichotomy result:', dichotomy(a, b, e, f))

    x1 = a + 0.381 * (b - a)
    x2 = a + 0.618 * (b - a)
    y1 = f(x1)
    y2 = f(x2)
    print('Golden Section result:', golden_section(a, b, x1, x2, y1, y2, e, f))

    print('Fibonacci result:', fibonacci(a, b, e, f))

    print('Parabola result:', parabola(a, 8.5 if f.__name__ == 'f2' else (a + b) / 2, b, e, f))

    print('Brent result:', brent(a, b, e, f))

    print()


# First Function
# for i in range(7):
#     e = 10 ** -i
a, b = -0.5, 0.5
e = 1e-5
test_func(a, b, e, f1)
# Second Function
a, b = 6, 9.9
e = 1e-5

test_func(a, b, e, f2)

# Third Function
a, b = 0, 2 * pi

e = 1e-5

test_func(a, b, e, f3)

# Fourth Function
a, b = 0, 1
e = 1e-5

test_func(a, b, e, f4)

# Fifth Function
a, b = 0.5, 2.5
e = 1e-5

test_func(a, b, e, f5)

a, b = -1, 2
e = 1e-5

test_func(a, b, e, f6)
