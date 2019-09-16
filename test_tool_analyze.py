import unittest
import numpy as np
import tool_analyze as ta


class TestToolAnalyze(unittest.TestCase):
    def test_stack_curves(self):
        x = []
        y = []
        x.append(np.array([1, 2, 3, 4, 5, 6, 7, 8]))
        y.append(np.ones(len(x[0]), dtype=np.float64))
        x.append(np.array([i for i in range(10)]))
        y.append(np.array([0, 1, 0, 0, 0, 0, 0, 2, 0, 0], dtype=np.float64))
        xnew, ynew = ta.stack_curves(x, y, 0.5, threshold=2)
        print("Y1:", y[0])
        print("Y2:", y[1])
        print("YS:", ynew)
        print("XS:", xnew)
        print(len(xnew) == len(ynew))

    def test_mass_bin(self):
        mbin = [0, 10, 20, 30, 40]
        R = [[1, 2, 3, 4, 5], [11, 12, 13, 14, 15]]
        Q = [[1, 1, 1, 1, 1], [2, 2, 2, 2, 2]]
        M = [40, 0]
        Rnew, Qnew = ta.mass_bin(R, Q, M, mbin)
        print('Rnew', Rnew)
        print('Qnew', Qnew)


if __name__ == "__main__":
    C = TestToolAnalyze()
    C.test_stack_curves()
    C.test_mass_bin()
    # unittest.main()
