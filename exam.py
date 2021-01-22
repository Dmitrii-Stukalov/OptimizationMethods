from math import sqrt


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

    def inverse(self):
        l, u = self.lu()
        b = [[1 if i == j else 0 for i in range(self.N)] for j in range(self.N)]
        inv = [[0 for _ in range(self.N)] for _ in range(self.N)]

        for i in range(self.N):
            y = self.forward(l, b[i])
            inv[i] = self.backward(u, y)

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


def test(A, max_k=10):
    x = [i + 1 for i in range(A.N)]
    for k in range(max_k):
        F = []
        for i in range(A.N):
            F.append(sum([A.to_matrix()[i][j] * x[j] for j in range(A.N)]))
        x_solve = A.lu_solve(F)

        A_noise = A.add_noise(k)
        F = []
        for i in range(A.N):
            F.append(sum([A_noise.to_matrix()[i][j] * x[j] for j in range(A_noise.N)]))
        x_noise = A.lu_solve(F)

        print('Expected:', x)
        print('Solved:', x_solve)
        print('With noise:', x_noise)
        print('Error',
              sqrt(sum([(x_solve[i] - x_noise[i]) ** 2 for i in range(A.N)]) / sum(map(lambda s: s ** 2, x_solve))))


A = [[9, 0, 0, 3, 1, 0, 1],
     [0, 11, 2, 1, 0, 0, 2],
     [0, 1, 10, 2, 0, 0, 0],
     [2, 1, 2, 9, 1, 0, 0],
     [1, 0, 0, 1, 12, 0, 1],
     [0, 0, 0, 0, 0, 8, 0],
     [2, 2, 0, 0, 3, 0, 8]]

B = CSRMatrix(A)

L, U = B.lu()
L.dot(U.to_matrix()).show(as_matrix=True)

B_inv = B.inverse()
B_inv.dot(B.to_matrix()).show(as_matrix=True)

A = [[1, 2, 3],
     [4, 5, 6],
     [7, 8, 0]]
b = [6, 9, -6]
analytic_ans = [-2, 1, 2]

A = CSRMatrix(A)
print(analytic_ans)
print(A.lu_solve(b))

A = [[1, 2, 3],
     [4, 5, 6],
     [7, 8, 0]]
A = CSRMatrix(A)
test(A)
