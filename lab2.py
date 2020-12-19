import operator
from itertools import cycle

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from prettytable import PrettyTable


def calc_func(func_index, x):
    if func_index == 1:
        return -5 * x ** 5 + 4 * x ** 4 - 12 * x ** 3 + 11 * x ** 2 - 2 * x + 1
    elif func_index == 2:
        return -3 * x * np.sin(0.75 * x) + np.exp(-2 * x)
    elif func_index == 3:
        return np.exp(3 * x) + 5 * np.exp(-2 * x)


def calc_der(func_index, x):
    if func_index == 1:
        return -25 * x ** 4 + 16 * x ** 3 - 36 * x ** 2 + 22 * x - 2
    elif func_index == 2:
        return -2 * np.exp(-2 * x) - 3 * np.sin(0.75 * x) - 2.25 * x * np.cos(0.75 * x)
    elif func_index == 3:
        return np.exp(-2 * x) * (3 * np.exp(5 * x) - 10)


def secant(x, w, f_x, f_w):
    return (x * f_w - w * f_x) / (f_w - f_x)


def brent(a, c, eps, func_index):
    x = w = v = (a + c) / 2
    f_x = f_w = f_v = calc_func(func_index, x)
    f_x_der = f_w_der = f_v_der = calc_der(func_index, x)
    d = e = c - a
    counter = 0
    prev_u = 0
    while True:
        counter += 1
        g, e = e, d
        u = None
        if x != w and f_x_der != f_w_der:
            u = secant(x, w, f_x_der, f_w_der)
            if a + eps <= u <= c - eps and abs(u - x) < g / 2:
                u = u
            else:
                u = None
        if x != v and f_x_der != f_v_der:
            u2 = secant(x, v, f_x_der, f_v_der)
            if a + eps <= u2 <= c - eps and abs(u2 - x) < g / 2:
                if u is not None and abs(u2 - x) < abs(u - x):
                    u = u2
        if u is None:
            if f_x_der > 0:
                u = (a + x) / 2
            else:
                u = (x + c) / 2
        if abs(u - x) < eps:
            u = x + np.sign(u - x) * eps
        d = abs(x - u)
        f_u = calc_func(func_index, u)
        f_u_der = calc_der(func_index, u)
        if f_u <= f_x:
            if u >= x:
                a = x
            else:
                c = x
            v, w, x = w, x, u
            f_v, f_w, f_x = f_w, f_x, f_u
            f_v_der, f_w_der, f_x_der = f_w_der, f_x_der, f_u_der
        else:
            if u >= x:
                c = u
            else:
                a = u
            if f_u <= f_w or w == x:
                v, w = w, u
                f_v, f_w = f_w, f_u
                f_v_der, f_w_der = f_w_der, f_u_der
            elif f_u <= f_v or v == x or v == w:
                v = u
                f_v = f_u
                f_v_der = f_u_der
        if counter > 1:
            if abs(prev_u - u) < eps:
                break
        prev_u = u
    return (a + c) / 2, counter


numbers = [1, 2, 3]
intervals = [[-0.5, 0.5], [0, np.pi * 2], [0, 1]]

table = PrettyTable()
table.field_names = ['function index', 'Min X', 'Min Y']

for number, interval in zip(numbers, intervals):
    x, counter = brent(interval[0], interval[1], 0.0001, number)
    table.add_row([number, x, '%.10f' % calc_func(number, x)])
print(table.get_string())

epsilons = np.linspace(1, 10, 10)
iterations = [brent(-0.5, 0.5, 10 ** -eps, 1)[1] for eps in epsilons]
plt.plot(np.log10(10 ** -epsilons), iterations)
plt.grid()
plt.ylabel('Количество итераций')
plt.xlabel('lg(eps)')
plt.title('Зависимость количества итераций от lg(eps), метод Брендта с производной')
plt.show()


def f1(x):
    return 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2


def f1_der(x, index):
    return {
        1: 2 * (200 * x[0] ** 3 - 200 * x[0] * x[1] + x[0] - 1),
        2: 200 * (x[1] - x[0] ** 2)
    }.get(index)


def f2(x):
    return (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2


def f2_der(x, index):
    return {
        1: 2 * (2 * x[0] ** 3 - 2 * x[0] * x[1] + x[0] - 1),
        2: 2 * (x[1] - x[0] ** 2)
    }.get(index)


def f3(x):
    return (1.5 - x[0] * (1 - x[1])) ** 2 + (2.25 - x[0] * (1 - x[1] ** 2)) ** 2 + (2.625 - x[0] * (1 - x[1] ** 3)) ** 2


def f3_der(x, index):
    return {
        1: 2 * x[0] * (
                x[1] ** 6 + x[1] ** 4 - 2 * x[1] ** 3 - x[1] ** 2 - 2 * x[1] + 3) + 5.25 * x[1] ** 3 + 4.5 * x[
               1] ** 2 + 3 * x[1] - 12.75,
        2: x[0] * (x[0] * (6 * x[1] ** 5 + 4 * x[1] ** 3 - 6 * x[1] ** 2 - 2 * x[1] - 2) + 15.75 * x[1] ** 2 + 9 * x[
            1] + 3)
    }.get(index)


def f4(x):
    return (x[0] + x[1]) ** 2 + 5 * (x[2] - x[3]) ** 2 + (x[1] - 2 * x[2]) ** 4 + 10 * (x[0] - x[3]) ** 4


def f4_der(x, index):
    return {
        1: 2 * (20 * (x[0] - x[3]) ** 3 + x[0] + x[1]),
        2: 2 * (x[0] + 2 * (x[1] - 2 * x[2]) ** 3 + x[1]),
        3: 10 * (x[2] - x[3]) - 8 * (x[1] - 2 * x[2]) ** 3,
        4: 10 * (-4 * (x[0] - x[3]) ** 3 + x[3] - x[2])
    }.get(index)


def get_func(index):
    return {1: f1,
            2: f2,
            3: f3,
            4: f4,
            }.get(index)


def get_func_der(index):
    return {1: f1_der,
            2: f2_der,
            3: f3_der,
            4: f4_der,
            }.get(index)


def get_func_dimension(index):
    return {1: 2,
            2: 2,
            3: 2,
            4: 4,
            }.get(index)


def get_bounds(index):
    return {1: (-1.0, 1.0),
            2: (-1.0, 1.0),
            3: (-3.0, 3.0),
            4: (-1.0, 1.0),
            }.get(index)


def init_approx(index):
    return {1: [-1.0, -1.0],
            2: [-1.0, -1.0],
            3: [-1.0, -1.0],
            4: [-1.0, -1.0, -1.0, -1.0],
            }.get(index)


def init_approx_ravine(index):
    return {1: ([1.5, 1.5], [1.45, 1.4]),
            2: ([0.5, 0.5], [0.45, 0.4]),
            3: ([0.5, 0.5], [0.45, 0.4]),
            4: ([0.5, 0.5, 0.5, 0.5], [0.45, 0.4, 0.5, 0.5]),
            }.get(index)


num_funcs = 4


def coordinate_descent(x, func_number, eps, learning_rate):
    coordinates = []
    n_arg = len(x)
    cycler = cycle([i + 1 for i in range(n_arg)])
    counter = 0
    while True:
        counter += 1
        index = next(cycler)
        d_x = get_func_der(func_number)(x, index)
        new_x = x.copy()
        new_x[index - 1] -= learning_rate * d_x
        coordinates.append([x, new_x])
        if np.linalg.norm(np.array(new_x) - np.array(x)) < eps:
            x = new_x
            break
        x = new_x
    return x, counter, coordinates


table = PrettyTable()
table.field_names = ['function index', 'Min X', 'Min Y']
for i in range(1, num_funcs + 1):
    x, counter, coordinates = coordinate_descent(init_approx(i), i, 0.00001, 0.002)
    table.add_row([i, x, '%.10f' % get_func(i)(x)])
print(table.get_string())

for i in range(1, num_funcs):
    delta = 0.025
    bounds = get_bounds(i)
    x = np.arange(bounds[0], bounds[1], delta)
    y = np.arange(bounds[0], bounds[1], delta)
    X, Y = np.meshgrid(x, y)
    Z = get_func(i)([X, Y])
    fig, ax = plt.subplots()
    cs = ax.contour(X, Y, Z)
    ax.clabel(cs)
    x, counter, coordinates = coordinate_descent(init_approx(i), i, 0.00001, 0.002)
    lc = LineCollection(coordinates)
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)
plt.show()

epsilons = np.linspace(1, 10, 10)
iterations = []
for eps in epsilons:
    x, counter, coordinates = coordinate_descent([-1.0, -1.0], 1, 10 ** -eps, 0.002)
    iterations.append(counter)
plt.plot(np.log10(10 ** -epsilons), iterations)
plt.grid()
plt.ylabel('Количество итераций')
plt.xlabel('lg(eps)')
plt.title('Зависимость количества итераци от lg(eps), метод покоординатного спуска')
plt.show()


def golden_section(a, b, func, eps):
    K = (1 + np.sqrt(5)) / 2
    x1 = b - (b - a) / K
    x2 = a + (b - a) / K
    f1 = func(x1)
    f2 = func(x2)
    counter = 0
    while abs(b - a) > eps:
        counter += 1
        if f1 < f2:
            b = x2
            f2, x2 = f1, x1
            x1 = b - (b - a) / K
            f1 = func(x1)
        else:
            a = x1
            f1, x1 = f2, x2
            x2 = a + (b - a) / K
            f2 = func(x2)
    return (a + b) / 2


def steepest_descent(x, func_index, eps):
    func = get_func(func_index)
    func_der = get_func_der(func_index)
    variable_number = get_func_dimension(func_index)
    x = [0] * variable_number
    points = [(x, func(x))]
    coordinates = []
    counter = 0
    while True:
        counter += 1
        grad = [func_der(x, i) for i in range(1, variable_number + 1)]
        new_x_func = lambda step: list(map(operator.sub, x, [step * grad[i] for i in range(len(grad))]))
        step_func = lambda step: func(new_x_func(step))
        result_step = golden_section(0, 10, step_func, 10 ** -5)
        new_x = new_x_func(result_step)
        coordinates.append([x, new_x])
        x = new_x
        points.append((x, func(x)))
        if abs(points[-1][1] - points[-2][1]) < eps:
            break
    return coordinates, points, counter


table = PrettyTable()
table.field_names = ['function index', 'Min X', 'Min Y']
for i in range(1, num_funcs + 1):
    coordinates, points, counter = steepest_descent(init_approx(i), i, 10 ** -9)
    table.add_row([i, points[-1][0], '%.15f' % points[-1][1]])
print(table.get_string())

for i in range(1, num_funcs):
    delta = 0.025
    bounds = get_bounds(i)
    x = np.arange(bounds[0], bounds[1], delta)
    y = np.arange(bounds[0], bounds[1], delta)
    X, Y = np.meshgrid(x, y)
    Z = get_func(i)([X, Y])
    fig, ax = plt.subplots()
    cs = ax.contour(X, Y, Z)
    ax.clabel(cs)
    coordinates, points, counter = steepest_descent(init_approx(i), i, 10 ** -9)
    lc = LineCollection(coordinates)
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)
plt.show()

iterations = []
for eps in epsilons:
    coordinates, points, counter = steepest_descent([0, 0], 1, 10 ** -eps)
    iterations.append(counter)
plt.plot(np.log10(10 ** -epsilons), iterations)
plt.grid()
plt.ylabel('Количество итераций')
plt.xlabel('lg(eps)')
plt.title('Зависимость количества итераций от lg(eps), метод наискорейшего спуска')


def projection(x, U, b):
    return b + U.dot(U.T).dot(x - b)


def steepest_descent_with_projection(x, func_index, eps):
    U = np.linalg.qr(np.random.normal(0, 1, (get_func_dimension(func_index), 100)))[0]
    b = np.random.normal(0, 1, get_func_dimension(func_index))
    func = get_func(func_index)
    func_der = get_func_der(func_index)
    variable_number = get_func_dimension(func_index)
    x = [0] * variable_number
    points = [(x, func(x))]
    coordinates = []
    counter = 0
    while True:
        counter += 1
        grad = [func_der(x, i) for i in range(1, variable_number + 1)]
        new_x_func = lambda step: list(map(operator.sub, x, [step * grad[i] for i in range(len(grad))]))
        step_func = lambda step: func(new_x_func(step))
        result_step = golden_section(0, 10, step_func, 10 ** -5)
        new_x = projection(np.array(new_x_func(result_step)), U, b)
        coordinates.append([x, new_x])
        x = new_x
        points.append((x, func(x)))
        if abs(points[-1][1] - points[-2][1]) < eps:
            break
    return coordinates, points, counter


table = PrettyTable()
table.field_names = ['function index', 'Min X', 'Min Y']
for i in range(1, num_funcs + 1):
    coordinates, points, counter = steepest_descent_with_projection(init_approx(i), i, 10 ** -9)
    table.add_row([i, points[-1][0], '%.15f' % points[-1][1]])
print(table.get_string())

for i in range(1, num_funcs):
    delta = 0.025
    bounds = get_bounds(i)
    x = np.arange(bounds[0], bounds[1], delta)
    y = np.arange(bounds[0], bounds[1], delta)
    X, Y = np.meshgrid(x, y)
    Z = get_func(i)([X, Y])
    fig, ax = plt.subplots()
    cs = ax.contour(X, Y, Z)
    ax.clabel(cs)
    coordinates, points, counter = steepest_descent_with_projection(init_approx(i), i, 10 ** -9)
    lc = LineCollection(coordinates)
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)
plt.show()

iterations = []
for eps in epsilons:
    coordinates, points, counter = steepest_descent([0, 0], 1, 10 ** -eps)
    iterations.append(counter)
plt.plot(np.log10(10 ** -epsilons), iterations)
plt.grid()
plt.ylabel('Количество итераций')
plt.xlabel('lg(eps)')
plt.title('Зависимость количества итераций от lg(eps), метод проекции градиента')
plt.show()


def step(v, func_number):
    x = []
    for i in range(len(v)):
        x.append(v[i] - 0.002 * get_func_der(func_number)(v, i + 1))
    return np.array(x)


def ravine_gradient(v_k, v_k_1, func_number, eps, C, h, max_iter):
    coordinates = []
    v_k = np.array(v_k)
    v_k_1 = np.array(v_k_1)
    x_k = step(v_k, func_number)
    x_k_1 = step(v_k_1, func_number)
    counter = 0
    while counter < max_iter:
        coordinates.append([x_k, x_k_1])
        counter += 1
        v_new = x_k_1 - (x_k_1 - x_k) / (np.linalg.norm(x_k_1 - x_k)) \
                * h * np.sign(np.array(get_func(func_number)(x_k_1) - np.array(get_func(func_number)(x_k))))
        x_new = step(v_new, func_number)
        cos_a = (np.dot(v_new - x_k_1, x_new - x_k_1)) / (
                np.linalg.norm(v_new - x_k_1) * (np.linalg.norm(x_new - x_k_1)))
        cos_a_1 = (np.dot(v_k_1 - x_k, x_k_1 - x_k)) / (np.linalg.norm(v_k_1 - x_k) * (np.linalg.norm(x_k_1 - x_k)))
        h *= C ** (cos_a - cos_a_1)
        v_k_1 = v_new
        x_k = x_k_1
        x_k_1 = x_new
        if np.linalg.norm(np.array(x_k_1) - np.array(x_k)) < eps:
            break
        if abs(get_func(func_number)(x_k_1) - get_func(func_number)(x_k)) < eps:
            break
    return x_k_1, counter, coordinates


table = PrettyTable()
table.field_names = ['function index', 'Min X', 'Min Y']
for i in range(1, num_funcs + 1):
    initialApprox = init_approx_ravine(i)
    x, counter, coordinates = ravine_gradient(initialApprox[0], initialApprox[1], i, 0.000001, 5, 0.0001, 5000)
    table.add_row([i, x, '%.15f' % get_func(i)(x)])
print(table.get_string())

initialApprox = []
for i in range(1, num_funcs):
    delta = 0.025
    bounds = get_bounds(i)
    x = np.arange(bounds[0], bounds[1], delta)
    y = np.arange(bounds[0], bounds[1], delta)
    X, Y = np.meshgrid(x, y)
    Z = get_func(i)([X, Y])
    fig, ax = plt.subplots()
    cs = ax.contour(X, Y, Z)
    ax.clabel(cs)
    initialApprox = init_approx_ravine(i)
    x, counter, coordinates = ravine_gradient(initialApprox[0], initialApprox[1], i, 0.00001, 3, 0.0001, 5000)
    lc = LineCollection(coordinates)
    ax.add_collection(lc)
    ax.autoscale()
    ax.margins(0.1)
plt.show()

iterations = []
for eps in epsilons:
    x, counter, coordinates = ravine_gradient(initialApprox[0], initialApprox[1], 1, 10 ** -eps, 10, 0.00002, 10000)
    iterations.append(counter)
plt.plot(np.log10(10 ** -epsilons), iterations)
plt.grid()
plt.ylabel('Количество итераций')
plt.xlabel('lg(eps)')
plt.title('Зависимость количества итераций от lg(eps), метод овражного градиента')
plt.show()

# def task - fix eps, diff init_approx
eps = 1e-5
initialApprox = init_approx_ravine(1)
x, counter, coordinates = ravine_gradient(initialApprox[0], initialApprox[1], 1, eps, 10, 0.00002, 10000)

print(counter)

