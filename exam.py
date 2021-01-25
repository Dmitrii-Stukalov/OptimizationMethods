from math import sqrt
from random import randint


class CSRMatrix(object):
    def __init__(self, matrix):
        self.N = len(matrix)
        self.data = []
        self.indptr = [1]
        self.indices = []

        for i in range(self.N):
            amount = 0
            for j in range(self.N):
                if matrix[i][j] != 0:
                    self.data.append(matrix[i][j])
                    self.indices.append(j)
                    amount += 1
            self.indptr.append(amount + self.indptr[-1])

    def to_matrix(self):
        matrix = [[0 for _ in range(self.N)] for _ in range(self.N)]
        data_ptr = 0
        col_ptr = 0
        for row in range(self.N):
            for j in range(self.indptr[row + 1] - self.indptr[row]):
                matrix[row][self.indices[col_ptr]] = self.data[data_ptr]
                data_ptr += 1
                col_ptr += 1

        return matrix

    def show(self, as_matrix=False, as_int=True):
        if as_matrix:
            for row in self.to_matrix():
                if as_int:
                    print([round(i) for i in row])
                else:
                    print(row)
        else:
            print('data:', self.data)
            print('indices:', self.indices)
            print('indptr:', self.indptr)

    def lu(self):
        l = [[0 for _ in range(self.N)] for _ in range(self.N)]
        u = [[0 for _ in range(self.N)] for _ in range(self.N)]

        for i in range(self.N):
            l[i][i] = 1

        a = self.to_matrix()
        for i in range(self.N):
            u[i][i] = a[i][i]
            for j in range(i + 1, self.N):
                l[j][i] = a[j][i] / u[i][i]
                u[i][j] = a[i][j]
            for j in range(i + 1, self.N):
                for k in range(i + 1, self.N):
                    a[j][k] -= l[j][i] * u[i][k]

        return CSRMatrix(l), CSRMatrix(u)

    def dot(self, B):
        B = B.to_matrix()
        answer = [[0 for _ in range(self.N)] for _ in range(self.N)]
        matrix = self.to_matrix()
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    answer[i][j] += matrix[i][k] * B[k][j]
        return CSRMatrix(answer)

    def forward(self, L, b):
        y = [0 for _ in range(self.N)]

        y[0] = b[0] / L.to_matrix()[0][0]
        for i in range(1, self.N):
            tmp = L.to_matrix()[i][:i]
            mul = [tmp[j] * y[j] for j in range(i)]
            y[i] = (b[i] - sum(mul)) / L.to_matrix()[i][i]
        return y

    def backward(self, U, y):
        x = [0 for _ in range(self.N)]

        x[-1] = y[-1] / U.to_matrix()[-1][-1]
        for i in range(self.N - 2, -1, -1):
            tmp = U.to_matrix()[i][i:]
            sub_x = x[i:]
            mul = [tmp[j] * sub_x[j] for j in range(len(tmp))]
            x[i] = (y[i] - sum(mul)) / U.to_matrix()[i][i]

        return x

    def get_column(self, column):
        return [row[column] for row in self.to_matrix()]

    def inverse(self):
        inv = zero(self.N).to_matrix()
        F = eye(self.N)
        for i in range(self.N):
            column = F.get_column(i)
            x = self.lu_solve(column)
            for j in range(self.N):
                inv[j][i] = x[j]

        return CSRMatrix(inv)

    def lu_solve(self, b):
        l, u = self.lu()
        y = self.forward(l, b)
        return self.backward(u, y)

    def add_noise(self, k):
        matrix = self.to_matrix()
        for i in range(self.N):
            matrix[i][i] += 10 ** -k
        return CSRMatrix(matrix)

    def sub(self, B):
        B = B.to_matrix()
        matrix = self.to_matrix()
        for i in range(self.N):
            for j in range(self.N):
                matrix[i][j] -= B[i][j]
        return CSRMatrix(matrix)


def norm(x):
    ans = 0
    for i in x:
        ans += i * i
    return sqrt(ans)


def test(A, max_k=15):
    A = A.to_matrix()
    for i in range(len(A)):
        A[i][i] = sum(A[i])
    A = CSRMatrix(A)
    x = [i + 1 for i in range(A.N)]
    for k in range(max_k):
        A_noise = A.add_noise(k)
        F = []
        for i in range(A.N):
            F.append(sum([A_noise.to_matrix()[i][j] * x[j] for j in range(A_noise.N)]))
        x_noise = A.lu_solve(F)

        print('Expected:', x)
        print('Solved:', x_noise)
        diff = [x[i] - x_noise[i] for i in range(len(x))]
        print('Error', norm(diff) / norm(x), 'k =', k)


def test2(max_k=15):
    for k in range(1, max_k):
        A = [[1 for _ in range(k)] for _ in range(k)]
        for i in range(k):
            for j in range(k):
                A[i][j] = 1 / (i + 1 + j + 1 - 1)
        A = CSRMatrix(A)
        x = [i + 1 for i in range(k)]
        F = []
        for i in range(A.N):
            F.append(sum([A.to_matrix()[i][j] * x[j] for j in range(A.N)]))
        x_solve = A.lu_solve(F)
        print('Expected:', x)
        print('Solved:', x_solve)
        diff = [x_solve[i] - x[i] for i in range(len(x_solve))]
        print('Error', norm(diff) / norm(x), 'k =', k)


def eye(k):
    matrix = [[0 for _ in range(k)] for _ in range(k)]
    for i in range(k):
        matrix[i][i] = 1
    return CSRMatrix(matrix)


def zero(k):
    return CSRMatrix([[0 for _ in range(k)] for _ in range(k)])


def block_lu(old_A, block_size):
    old_A = old_A.to_matrix()
    n = len(old_A)
    X = [[0 for _ in range(block_size)] for _ in range(block_size)]
    Y = [[0 for _ in range(n - block_size)] for _ in range(block_size)]
    Z = [[0 for _ in range(block_size)] for _ in range(n - block_size)]
    W = [[0 for _ in range(n - block_size)] for _ in range(n - block_size)]
    for i in range(n):
        for j in range(n):
            if i < block_size:
                if j < block_size:
                    X[i][j] = old_A[i][j]
                else:
                    Y[i][j - block_size] = old_A[i][j]
            else:
                if j < block_size:
                    Z[i - block_size][j] = old_A[i][j]
                else:
                    W[i - block_size][j - block_size] = old_A[i][j]
    X = CSRMatrix(X)
    Y = CSRMatrix(Y)
    Z = CSRMatrix(Z)
    W = CSRMatrix(W)

    L11, U11 = X.lu()
    U12 = L11.inverse().dot(Y)
    L12 = zero(Y.N)
    L21 = Z.dot(U11.inverse())
    U21 = zero(Z.N)

    L22, U22 = W.sub(L21.dot(U12)).lu()

    L = [[0 for _ in range(n)] for _ in range(n)]
    U = [[0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(n):
            if i < block_size:
                if j < block_size:
                    L[i][j] = L11.to_matrix()[i][j]
                    U[i][j] = U11.to_matrix()[i][j]
                else:
                    L[i][j] = L12.to_matrix()[i][j - block_size]
                    U[i][j] = U12.to_matrix()[i][j - block_size]
            else:
                if j < block_size:
                    L[i][j] = L21.to_matrix()[i - block_size][j]
                    U[i][j] = U21.to_matrix()[i - block_size][j]
                else:
                    L[i][j] = L22.to_matrix()[i - block_size][j - block_size]
                    U[i][j] = U22.to_matrix()[i - block_size][j - block_size]

    return CSRMatrix(L), CSRMatrix(U)


def generate_test_matrix(k):
    return CSRMatrix([[randint(-10, 10) for _ in range(k)] for _ in range(k)])


def generate_special_test_matrix(k):
    return CSRMatrix([[randint(-4, 0) for _ in range(k)] for _ in range(k)])


A = [[9, 0, 0, 3, 1, 0, 1],
     [0, 11, 2, 1, 0, 0, 2],
     [0, 1, 10, 2, 0, 0, 0],
     [2, 1, 2, 9, 1, 0, 0],
     [1, 0, 0, 1, 12, 0, 1],
     [0, 0, 0, 0, 0, 8, 0],
     [2, 2, 0, 0, 3, 0, 8]]

A = CSRMatrix(A)

L, U = A.lu()
L.dot(U).show(as_matrix=True)
print()

A_inv = A.inverse()
A_inv.dot(A).show(as_matrix=True)
print()

A = [[1, 2, 3],
     [4, 5, 6],
     [7, 8, 0]]
b = [6, 9, -6]
analytic_ans = [-2, 1, 2]

A = CSRMatrix(A)
print('Analytical answer:', analytic_ans)
print('Solved:', A.lu_solve(b))
print()

A = generate_special_test_matrix(4)
test(A)
test2()

for k in range(1, 8):
    A = generate_test_matrix(2 * k)
    A.show(as_matrix=True)
    print()
    block_size = randint(1, 2 * k)
    print('Block size:', block_size)
    L, U = block_lu(A, k)
    L.show(as_matrix=True)
    print()
    U.show(as_matrix=True)
    print()
    L.dot(U).show(as_matrix=True)
    print()
