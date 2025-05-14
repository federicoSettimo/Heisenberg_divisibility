#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
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
static Eigen::Matrix2cd sigma_p {{0,1},{0,0}};
static Eigen::Matrix2cd sigma_m {{0,0},{1,0}};
static Eigen::Matrix2cd id {{1,0},{0,1}};

static Eigen::Vector2cd ground_state {{0.,1.}};
static Eigen::Vector2cd excited_state {{1.,0.}};
static Eigen::Vector2cd plus_state {{1./sqrt(2.),1./sqrt(2.)}};
static Eigen::Vector2cd minus_state {{1./sqrt(2.),-1./sqrt(2.)}};

// Params
double tmax_rates = 2.2, dt = 0.001, tmax_norm = 5.;

double lambda_T (double t) {return .5*sin(t);}
double lambda_z (double t) {return exp(-t);}
double lambda (double t) {return exp(-2*t);}

double gamma_p (double t) {return 0.5 + 0.25*cos(t) + 0.25*sin(t);}
double gamma_m (double t) {return 0.5 - 0.25*cos(t) - 0.25*sin(t);}
double gamma_z (double t) {return 0.;}

double xi_p (double t) {return .5 + .25*exp(t)*cos(t);}
double xi_m (double t) {return .5 - .25*exp(t)*cos(t);}
double xi_z (double t) {return 0.;}

Matrix2cd PhiDag (const Matrix2cd& X, double t) {
    return .5*(
        ((id + lambda_T(t)*sigma_z)*X).trace() * id +
        lambda(t)*(X*sigma_x).trace() * sigma_x +
        lambda(t)*(X*sigma_y).trace() * sigma_y +
        lambda_z(t)*(X*sigma_z).trace() * sigma_z
    );
}

Matrix2cd proj (const Vector2cd& v) {return v*v.adjoint();}

double norm_infty (const Matrix2cd &X) {
    ComplexEigenSolver<Matrix2cd> eigs;
    eigs.compute(X);
    double l1 = abs(eigs.eigenvalues()(0)), l2 = abs(eigs.eigenvalues()(1));
    return l1 > l2 ? l1 : l2;
}

int main () {
    ofstream out_rates("rates.dat"), out_norm("norm.dat");
    out_rates << tmax_rates << endl << dt << endl;
    out_norm << tmax_norm << endl << dt << endl;

    // Plotting the rates
    for (double t = 0.; t < tmax_rates; t += dt) {
        out_rates << gamma_p(t) << " " << gamma_m(t) << " " << gamma_z(t) << endl;
        out_rates << xi_p(t) << " " << xi_m(t) << " " << xi_z(t) << endl;
    }

    // Now the norm
    // Initial effects - projector on pure states
    Matrix2cd E = proj(ground_state), F = proj(excited_state);
    for (double t = 0.; t < tmax_norm; t += dt) {
        Matrix2cd Et = PhiDag(E,t), Ft = PhiDag(F,t);
        out_norm << norm_infty(Et) << " " << norm_infty(Ft) << " " << norm_infty(Et-Ft) << endl;
    }

    return 0;
}