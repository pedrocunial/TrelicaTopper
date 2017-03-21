import numpy as np


class NumericMethods:

    def gauss_seidel(ite, tol, K, F):
        n = len(K)
        x = np.zeros(n)
        y = np.zeros(n)
        i = 0
        err = tol + 1
        while i < ite and err > tol:
            err = 0
            for j in range(n):
                y[j] = x[j]
                sum = 0
                for z in range(n):
                    if j != z:
                        sum += K[j][z] * x[z]
                x[j] = (F[j] - sum) / K[j][j]
                _err = abs((x[j] - y[j]) / x[j]) * 100
                if _err > err:
                    err = _err
            i += 1
        return x, err, i

    def jacobi(ite, tol, K, F):
        n = len(K)
        x = np.zeros(n)
        y = np.zeros(n)
        tmp = np.zeros(n)
        i = 0
        err = tol + 1
        while i < ite and err > tol:
            err = 0
            for k in range(n):
                tmp[k] = x[k]
            for j in range(n):
                y[j] = x[j]
                sum = 0
                for z in range(n):
                    if j != z:
                        sum += K[j][z] * tmp[z]
                x[j] = (F[j] - sum) / K[j][j]
                _err = abs((x[j] - y[j]) / x[j]) * 100
                if _err > err:
                    err = _err
            i += 1
        return x, err, i



if __name__ == "__main__":
    K = [
        [10, 2, -3, 2],
        [2, -15, 3, -2],
        [1, -3, 20, 2],
        [2, 2, -1, 30]
    ]
    F = [
        32,
        -59,
        -38,
        160
    ]

    x, err, ite = NumericMethods.gauss_seidel(4, 0, K, F)
    print("Gauss-Seidel 4 iterações")
    print(x)
    print(err)
    print(ite)

    x, err, ite = NumericMethods.jacobi(4, 0, K, F)
    print("Jacobi 4 iterações")
    print(x)
    print(err)
    print(ite)

    x, err, ite = NumericMethods.gauss_seidel(20, 0, K, F)
    print("Gauss-Seidel 20 iterações")
    print(x)
    print(err)
    print(ite)

    x, err, ite = NumericMethods.jacobi(20, 0, K, F)
    print("Jacobi 20 iterações")
    print(x)
    print(err)
    print(ite)
