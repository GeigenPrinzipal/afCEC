#ifndef UTILS
#define UTILS

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma ;

//   -*-   -*-   -*-

bool CholeskyRankOneUpdate(mat &L, vec v) {
    double *vPtr = v.memptr();
    for (int i = 0; i < L.n_rows; ++i) {
        double *LPtr = L.colptr(i);
        double r = sqrt((LPtr[i] * LPtr[i]) + (vPtr[i] * vPtr[i]));
        double c = LPtr[i] / r;
        double s = -vPtr[i] / r;
        LPtr[i] = r;
        for (int j = i + 1; j < L.n_cols; ++j) {
            double tmp1 = (LPtr[j] * c) - (vPtr[j] * s);
            double tmp2 = (LPtr[j] * s) + (vPtr[j] * c);
            LPtr[j] = tmp1;
            vPtr[j] = tmp2;
        }
    }
    return true;
}

//   -*-   -*-   -*-

bool CholeskyRankOneDowndate(mat &L, vec v) {
    double *vPtr = v.memptr();
    for (int i = 0; i < L.n_rows; ++i) {
        double *LPtr = L.colptr(i);
        double r = sqrt((LPtr[i] * LPtr[i]) - (vPtr[i] * vPtr[i]));
        double c = LPtr[i] / r;
        double s = -vPtr[i] / r;
        LPtr[i] = r;
        for (int j = i + 1; j < L.n_cols; ++j) {
            double tmp1 = (LPtr[j] * c) + (vPtr[j] * s);
            double tmp2 = (LPtr[j] * s) + (vPtr[j] * c);
            LPtr[j] = tmp1;
            vPtr[j] = tmp2;
        }
    }
    return true;
}

//   -*-   -*-   -*-

#endif
