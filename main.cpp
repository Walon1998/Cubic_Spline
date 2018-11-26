// to compile run: g++ -std=gnu++11 cubic_spline.cpp -lmgl
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <mgl2/mgl.h>

using namespace Eigen;
using namespace std;

MatrixXd cubicSpline(const VectorXd &T, const VectorXd &Y) {
    // returns the matrix representing the spline interpolating the data
    // with abscissae T and ordinatae Y. Each column represents the coefficients
    // of the cubic polynomial on a subinterval.
    // Assumes T is sorted, has no repeated elements and T.size() == Y.size().

    long n = T.size() - 1; // T and Y have length n+1

    VectorXd h(n);

    for (int i = 0; i < n; ++i) {
        h(i) = T(i + 1) - T(i);
    }

    MatrixXd A = MatrixXd::Zero(n - 1, n - 1);

    for (int j = 0; j < n - 1; ++j) {
        A(j, j) = (h(j) + h(j + 1)) / 3;
    }

    for (int j = 0; j < n - 2; ++j) {
        A(j + 1, j) = h(j + 1) / 6;
    }

    for (int j = 0; j < n - 2; ++j) {
        A(j, j + 1) = h(j + 1) / 6;
    }


    VectorXd rh(n);

    for (int i = 0; i < n; ++i) {
        rh(i) = (Y(i + 1) - Y(i)) / h(i);
    }

    VectorXd r(n - 1);

    for (int k = 0; k < n - 1; ++k) {
        r(k) = rh(k + 1) - rh(k);
    }

    VectorXd result(n - 1);
    result = A.partialPivLu().solve(r);

    VectorXd sigma(n + 1);
    sigma(0) = 0;
    sigma(n) = 0;
    for (int l = 1; l < n; ++l) {
        sigma(l) = result(l - 1);
    }

    VectorXd sigmah(n);

    for (int i = 0; i < n; ++i) {
        sigmah(i) = sigma(i + 1) + 2 * sigma(i);
    }

    VectorXd sigmah2(n);

    for (int i = 0; i < n; ++i) {
        sigmah2(i) = sigma(i + 1) - sigma(i);
    }

    MatrixXd spline(4, n);


    for (int m = 0; m < n; ++m) {
        spline(0, m) = Y(m);
        spline(1, m) = rh(m) - h(m) * ((sigmah(m)) / 6);
        spline(2, m) = sigma(m) / 2;
        spline(3, m) = sigmah2(m) / (6 * h(m));
    }


    return spline;
}

VectorXd evalCubicSpline(const MatrixXd &S, const VectorXd &T, const VectorXd &evalT) {
    // Returns the values of the spline S calculated in the points X.
    // Assumes T is sorted, with no repetetions.

    long n = evalT.size();

    VectorXd out(n);
    long m = T.size();


    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m - 1; j++) {
            if (evalT(i) < T(j + 1) || j == m - 2) {
                double value = evalT(i) - T(j);
                double a = value * S(1, j);
                double b = value * value * S(2, j);
                double c = value * value * value * S(3, j);

                double ab = a + b;
                double cd = c + S(0, j);
                out(i) = ab + cd;


                break;
            }

        }
    }


    cout << out << endl;

    return out;
}

int main() {
    // tests
    VectorXd T(9);
    VectorXd Y(9);
    T << 0, 0.4802, 0.7634, 1, 1.232, 1.407, 1.585, 1.879, 2;
    Y << 0., 0.338, 0.7456, 0, -1.234, 0, 1.62, -2.123, 0;

    int len = 1 << 9;
    VectorXd evalT = VectorXd::LinSpaced(len, T(0), T(T.size() - 1));

    VectorXd evalSpline = evalCubicSpline(cubicSpline(T, Y), T, evalT);

    mglData datx, daty;
    datx.Link(evalT.data(), len);
    daty.Link(evalSpline.data(), len);
    mglGraph gr;
    gr.SetRanges(0, 2, -3, 3);
    gr.Plot(datx, daty, "0");
    gr.WriteFrame("spline.eps");
}
	
