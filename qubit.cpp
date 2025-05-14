#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/MatrixFunctions>
#include <complex>
#include <vector>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace Eigen;

static complex<double> I(0,1), one(1,0), zero(0,0);
static Eigen::Matrix2cd sigma_x {{0,1},{1,0}};
static Eigen::Matrix2cd sigma_y {{0,-I},{I,0}};
static Eigen::Matrix2cd sigma_z {{1,0},{0,-1}};
static Eigen::Matrix2cd id {{1,0},{0,1}};

// Parameters
double tmax = 3., dt = .01, dlambda = .0005;

Matrix2cd BlochToMatrix (const Vector3d &x) {
    return .5*(id + x(0)*sigma_x + x(1)*sigma_y + x(2)*sigma_z);
}

Vector4d MatrixToBloch (const Matrix2cd &X) {
    Vector4d x {{real(X.trace()), real((X*sigma_x).trace()), real((X*sigma_y).trace()), real((X*sigma_z).trace())}};
    return x;
}

// Dynamics
double l1 (double t) {
    double lf = .4;
    return t <1. ? (1.-lf)*(1.-t)*(1.-t) + lf : lf;
}
double l2 (double t) {return l1(t);}
double l3 (double t) {return 1.;}

double a (double t) {return 0.;}
double b (double t) {return 0.;}
double c (double t) {return t < 1. ? 0. : (t-1.);}

Matrix3d D (double t) {return DiagonalMatrix<double, 3>(l1(t), l2(t), l3(t));}
Matrix3d A (double t) {
    Matrix3d res;
    res << 0., a(t), b(t), -a(t), 0., c(t), -b(t), -c(t), 0.;
    return res;
}
Matrix3d O (double t) {return A(t).exp();}
Matrix3d Lambda (double t) {return D(t)*O(t).transpose();}
Matrix3d Lambda_S (double t) {return O(t)*D(t);}


// Distances
double norm_infty (const Matrix2cd &X) {
    ComplexEigenSolver<Matrix2cd> eigs;
    eigs.compute(X);
    double l1 = abs(eigs.eigenvalues()(0)), l2 = abs(eigs.eigenvalues()(1));
    return l1 > l2 ? l1 : l2;
}

double OD (const Vector3d &r, const Vector3d &s) {
    return norm_infty(BlochToMatrix(r) - BlochToMatrix(s));
}

double TD (const Vector3d &r, const Vector3d &s) {
    return .5*(r-s).norm();
}

// Sharpness
double sharpness (const Matrix2cd &X) {
    return norm_infty(X) + norm_infty(id - X) - 1.;
}

// Incompatibility
double MinkProd (const Vector4d& x, const Vector4d& y) {
    return x(0)*y(0) - x(1)*y(1) - x(2)*y(2) - x(3)*y(3);
}

double MinkSquare (const Vector4d& x) {return MinkProd(x,x);}

Vector4d perpPOVM (const Vector4d& x) {
    Vector4d xperp {{2-x(0), -x(1), -x(2), -x(3)}};
    return xperp;
}

double comp (const Vector4d& x, const Vector4d& y) {
    Vector4d xperp = perpPOVM(x), yperp = perpPOVM(y);
    return sqrt(MinkSquare(x)*MinkSquare(xperp)*MinkSquare(y)*MinkSquare(yperp)) - 
        MinkProd(x,xperp)*MinkProd(y,yperp) +
        MinkProd(x,yperp)*MinkProd(xperp,y) +
        MinkProd(x,y)*MinkProd(xperp,yperp);
}

Vector4d N (const Vector4d &x, double l, double b = 0.) {
    Vector4d Nx {{(1.-l)*x(0) + l*(1.+b), (1.-l)*x(1), (1.-l)*x(2), (1.-l)*x(3)}};
    return Nx;
}

double Inc (const Vector3d &x, const Vector3d &y, double b = 0.) {
    Vector4d xprime {{1., x(0), x(1), x(2)}}, yprime {{1., y(0), y(1), y(2)}};
    for (double l = 0.; l <= .5; l += dlambda) {
        if (comp(N(xprime,l,b), N(yprime,l,b)) >= 0.)
            return l;
    }
    return .5;
}

double Inc_steering (const Vector3d &x, const Vector3d &y) {
    //Vector4d xprime {{1., x(0), x(1), x(2)}}, yprime {{1., y(0), y(1), y(2)}};
    for (double l = 0.; l <= 1.; l += dlambda) {
        Matrix2cd X = BlochToMatrix(x), Y = BlochToMatrix(y);
        if (comp(MatrixToBloch((1.-l)*X + .5*X.trace()*l*id), MatrixToBloch((1.-l)*Y + .5*Y.trace()*l*id)) >= 0.)
            return l;
    }
    return 1.;
}

double rand01 () {return rand()/((double)RAND_MAX);}
double min2 (double x, double y) {return x < y ? x : y;}

int main () {
    Vector3d plus_x {{1., 0., 0.}}, plus_z {{0., 0., 1.}}, plus_y {{0., 1., 0.}},
        minus_x {{-1., 0., 0.}}, minus_z {{0., 0., -1.}}, minus_y {{0., -1., 0.}};

    ofstream params("params.txt"), out_norms("norms.txt"), out_s_i("sh_inc.txt"), out_funcs("functions.txt");

    params << tmax << endl << dt << endl;

    for (double t = 0.; t < tmax; t += dt) {
        Vector3d m_y_t = Lambda(t)*minus_y, p_y_t = Lambda(t)*plus_y, p_x_t = Lambda(t)*plus_x;
        Vector3d m_y_t_S = Lambda_S(t)*minus_y, p_y_t_S = Lambda_S(t)*plus_y;

        out_norms << OD(m_y_t, p_y_t) << " " << TD(m_y_t_S, p_y_t_S) << endl;
        out_s_i << sharpness(BlochToMatrix(p_y_t)) << " " << Inc(p_y_t, p_x_t, 0.) << " " << Inc_steering(p_y_t, p_x_t) << endl;
        out_funcs << l1(t) << " " << c(t) << endl;
    }

    return 0;
}