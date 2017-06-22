#ifndef C_AFCEC_HARTIGAN
#define C_AFCEC_HARTIGAN

#include <stdio.h>
#include <conio.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include ".\Utils.h"

using namespace arma;

//   -*-   -*-   -*-

class CafCECHartigan {
    public:
        CafCECHartigan(
            const mat &points,
            int maxClusters,
            ivec labels,
            double cardMin,
            int minIterations,
            int maxIterations,
            double dEMin,
            const mat &values,
            bool interactive
        ) :
            points(points),
            maxClusters(maxClusters),
            labels(labels),
            cardMin(cardMin),
            minIterations(minIterations),
            maxIterations(maxIterations),
            dEMin(dEMin),
            values(values),
            interactive(interactive)
        {
            Initialize();
            SClusterData bestDataExcl(*this);
            SClusterData data(*this);
            SClusterData bestDataIncl(*this);
            numberOfIterations = 0;
            EOld = INFINITY;
            ENew = 0.0;
            for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) ENew += cls[*it].E;
            printf("%d, %lf\n", activeClusters.size(), ENew); // !!!
            while ((numberOfIterations < minIterations) || ((numberOfIterations < maxIterations) && (ENew - EOld < -dEMin))) {
                if (interactive) UpdateResult();
                for (int i = 0; i < pointsNum; ++i) {
                    int cl = labels[i];
                    if (activeClusters.find(cl) != activeClusters.end()) {
                        double dEMinExcl;
                        try {
                            bestDataExcl = cls[cl];
                            bestDataExcl.Downdate(i);
                            dEMinExcl = bestDataExcl.E - cls[cl].E;
                        } catch (...) {
                            RemoveCluster(cl);
                            activeClusters.erase(activeClusters.find(cl));
                        }
                        if (activeClusters.find(cl) != activeClusters.end()) {
                            double dEMinIncl = INFINITY;
                            int bestClIncl;
                            for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
                                if (*it != cl) {
                                    data = cls[*it];
                                    data.Update(i);
                                    if (data.E - cls[*it].E < dEMinIncl) {
                                        dEMinIncl = data.E - cls[*it].E;
                                        bestClIncl = *it;
                                        bestDataIncl = data;
                                    }
                                }
                            }
                            if (dEMinExcl + dEMinIncl < 0.0) {
                                cls[cl] = bestDataExcl;
                                cls[bestClIncl] = bestDataIncl;
                                labels[i] = bestClIncl;
                                if ((cls[cl].card <= dim) || (((double)cls[cl].card) / pointsNum < cardMin)) {
                                    RemoveCluster(cl);
                                    activeClusters.erase(activeClusters.find(cl));
                                }
                            }
                        }
                    }
                }
                EOld = ENew;
                for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) ENew += cls[*it].E;
                ++numberOfIterations;
            }
            UpdateResult();
        }

        ~CafCECHartigan() {
            delete[] sigma;
        }
    public:
        struct SClusterData {
            const CafCECHartigan &afCECHartigan;

            int card;
            double prob;
            vec m;
            mat *sigmaL;
            mat *A;
            mat *AL;
            vec *b;
            double *sumOfSquares;
            int axis;
            vec coeffs;
            double var;
            double E;

            //   -*-   -*-   -*-

            SClusterData(const CafCECHartigan &afCECHartigan) :
                afCECHartigan(afCECHartigan),

                card(0),
                prob(0.0),
                axis(0),
                var(0.0),
                E(0.0)
            {
                m.zeros(afCECHartigan.dim);
                sigmaL = new mat[afCECHartigan.dim];
                A = new mat[afCECHartigan.dim];
                AL = new mat[afCECHartigan.dim];
                b = new vec[afCECHartigan.dim];
                sumOfSquares = new double[afCECHartigan.dim];
                for (int i = 0; i < afCECHartigan.dim; ++i) {
                    sigmaL[i].zeros(afCECHartigan.dim, afCECHartigan.dim);
                    A[i].zeros(afCECHartigan.dimA, afCECHartigan.dimA);
                    AL[i].zeros(afCECHartigan.dimA, afCECHartigan.dimA);
                    b[i].zeros(afCECHartigan.dimA);
                    sumOfSquares[i] = 0.0f;
                }
                coeffs.zeros(afCECHartigan.dimA);
            }

            //   -*-   -*-   -*-

            SClusterData(const SClusterData &data) :
                afCECHartigan(data.afCECHartigan),

                card(data.card),
                prob(data.prob),
                axis(data.axis),
                var(data.var),
                E(data.E)
            {
                /*sigmaL(NULL),
                A(NULL),
                AL(NULL),
                b(NULL),
                sumOfSquares(NULL),*/

                m = data.m;
                sigmaL = new mat[afCECHartigan.dim];
                A = new mat[afCECHartigan.dim];
                AL = new mat[afCECHartigan.dim];
                b = new vec[afCECHartigan.dim];
                sumOfSquares = new double[afCECHartigan.dim];
                for (int i = 0; i < afCECHartigan.dim; ++i) {
                    sigmaL[i] = data.sigmaL[i];
                    A[i] = data.A[i];
                    AL[i] = data.AL[i];
                    b[i] = data.b[i];
                    sumOfSquares[i] = data.sumOfSquares[i];
                }
                coeffs = data.coeffs;
            }

            //   -*-   -*-   -*-

            SClusterData operator=(const SClusterData &data) {
                card = data.card;
                prob = data.prob;
                m = data.m;
                for (int i = 0; i < afCECHartigan.dim; ++i) {
                    sigmaL[i] = data.sigmaL[i];
                    A[i] = data.A[i];
                    AL[i] = data.AL[i];
                    b[i] = data.b[i];
                    sumOfSquares[i] = data.sumOfSquares[i];
                }
                axis = data.axis;
                coeffs = data.coeffs;
                var = data.var;
                E = data.E;
            }

            //   -*-   -*-   -*-

            ~SClusterData() {
                for (int i = 0; i < afCECHartigan.dim; ++i) {
                    if (sigmaL != NULL) delete[] sigmaL;
                    if (A != NULL) delete[] A;
                    if (AL != NULL) delete[] AL;
                    if (b != NULL) delete[] b;
                    if (sumOfSquares != NULL) delete sumOfSquares;
                }
            }

            //   -*-   -*-   -*-

            void Update(int point) {
                double dEMin = INFINITY;
                int bestAxis;
                vec bestCoeffs;
                double bestVar;
                double EMin;

                ++card;
                prob = ((double)card) / afCECHartigan.pointsNum;
                vec v = (afCECHartigan.points.submat(1, point, size(afCECHartigan.dim - 1, 1)) - m.subvec(1, afCECHartigan.dim - 1)) / sqrt((double)card);
                for (int i = 0; i < afCECHartigan.dim; ++i) {
                    CholeskyRankOneUpdate(sigmaL[i], v);
                    sigmaL[i] *= sqrt((card - 1.0) / card);
                    double det = 1.0;
                    for (int j = 0; j < afCECHartigan.dim - 1; ++j) det *= sigmaL[i].at(j, j);

                    vec vA = afCECHartigan.values.submat(i * afCECHartigan.dimA, point, size(afCECHartigan.dimA, 1));
                    CholeskyRankOneUpdate(AL[i], vA);
                    b[i] = b[i] + (afCECHartigan.points.at(i, point) * vA);
                    vec y = solve(trimatl(AL[i]), b[i]);
                    vec x = solve(trimatu(AL[i].t()), y);
                    A[i] = A[i] + (vA * vA.t());
                    sumOfSquares[i] = sumOfSquares[i] + (afCECHartigan.points.at(i, point) * afCECHartigan.points.at(i, point));
                    double var = sumOfSquares[i];
                    var -= ((vec)(2.0 * x.t() * b[i]))[0];
                    var += ((mat)(x.t() * A[i] * x)).at(0, 0);
                    var /= card;

                    double ETmp = prob * (-log(prob) + (0.5 * ((afCECHartigan.dim * log(2.0 * M_PI * M_E)) + log(det * det * var))));
                    if (!is_finite(ETmp)) throw "A numerical overflow occurred.";
                    if (ETmp - E < dEMin) {
                        dEMin = ETmp - E;
                        bestAxis = i; bestCoeffs = x; bestVar = var; EMin = ETmp;
                    }

                    if (i < afCECHartigan.dim - 1) v[i] = (afCECHartigan.points.at(i, point) - m[i]) / sqrt((double)card);
                }
                m = ((m * (card - 1)) + afCECHartigan.points.col(point)) / card;
                axis = bestAxis;
                coeffs = bestCoeffs;
                var = bestVar;
                E = EMin;
            }

            //   -*-   -*-   -*-

            SClusterData Downdate(int point) {
                double dEMin = INFINITY;
                int bestAxis;
                vec bestCoeffs;
                double bestVar;
                double EMin;

                --card;
                prob = ((double)card) / afCECHartigan.pointsNum;
                vec v = (afCECHartigan.points.submat(1, point, size(afCECHartigan.dim - 1, 1)) - m.subvec(1, afCECHartigan.dim - 1)) / sqrt((double)card);
                for (int i = 0; i < afCECHartigan.dim; ++i) {
                    CholeskyRankOneDowndate(sigmaL[i], v);
                    sigmaL[i] *= sqrt((card + 1.0) / card);
                    double det = 1.0;
                    for (int j = 0; j < afCECHartigan.dim - 1; ++j) det *= sigmaL[i].at(j, j);

                    vec vA = afCECHartigan.values.submat(i * afCECHartigan.dimA, point, size(afCECHartigan.dimA, 1));
                    CholeskyRankOneDowndate(AL[i], vA);
                    b[i] = b[i] - (afCECHartigan.points.at(i, point) * vA);
                    vec y = solve(trimatl(AL[i]), b[i]);
                    vec x = solve(trimatu(AL[i].t()), y);
                    A[i] = A[i] - (vA * vA.t());
                    sumOfSquares[i] = sumOfSquares[i] - (afCECHartigan.points.at(i, point) * afCECHartigan.points.at(i, point));
                    double var = sumOfSquares[i];
                    var -= ((vec)(2.0 * x.t() * b[i]))[0];
                    var += ((mat)(x.t() * A[i] * x)).at(0, 0);
                    var /= card;

                    double ETmp = prob * (-log(prob) + (0.5 * ((afCECHartigan.dim * log(2.0 * M_PI * M_E)) + log(det * det * var))));
                    if (!is_finite(ETmp)) throw "A numerical overflow occurred.";
                    if (ETmp - E < dEMin) {
                        dEMin = ETmp - E;
                        bestAxis = i; bestCoeffs = x; bestVar = var; EMin = ETmp;
                    }

                    if (i < afCECHartigan.dim - 1) v[i] = (afCECHartigan.points.at(i, point) - m[i]) / sqrt((double)card);
                }
                m = ((m * (card + 1)) - afCECHartigan.points.col(point)) / card;
                axis = bestAxis;
                coeffs = bestCoeffs;
                var = bestVar;
                E = EMin;
            }
        };

        //   -*-   -*-   -*-

        const mat &points;
        int maxClusters;
        ivec labels;
        double cardMin;
        int minIterations;
        int maxIterations;
        double dEMin;
        const mat &values;
        bool interactive;

        Rcpp::List res;
        int pointsNum;
        int dim;
        int dimA;
        std::vector<SClusterData> cls;
        std::set<int> activeClusters;
        std::set<int> inactiveClusters;
        mat *sigma;
        int numberOfIterations;
        double EOld;
        double ENew;

        void Initialize() {
            pointsNum = points.n_cols;
            dim = points.n_rows;
            dimA = values.n_rows / dim;

            SClusterData tmp(*this);
            cls = std::vector<SClusterData>(maxClusters, tmp);

            for (int i = 0; i < pointsNum; ++i) ++cls[labels[i]].card;
            for (int i = 0; i < maxClusters; ++i) {
                if ((cls[i].card >= dim + 1) && (((double)cls[i].card) / pointsNum >= cardMin)) {
                    cls[i].prob = ((double)cls[i].card) / pointsNum;
                    activeClusters.insert(activeClusters.end(), i);
                } else {
                    if (cls[i].card > 0) inactiveClusters.insert(inactiveClusters.end(), i);
                }
            }

            sigma = new mat[maxClusters];
            for (int i = 0; i < pointsNum; ++i) {
                int cl = labels[i];
                if (activeClusters.find(cl) != activeClusters.end()) {
                    vec v = points.col(i);
                    cls[cl].m += v;
                    for (int j = 0; j < dim; ++j) {
                        vec vA = values.submat(j * dimA, i, size(dimA, 1));
                        cls[cl].A[j] += vA * vA.t();
                        cls[cl].b[j] += v[j] * vA;
                        cls[cl].sumOfSquares[j] += v[j] * v[j];
                    }
                }
            }
            for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
                cls[*it].m /= cls[*it].card;
                sigma[*it].zeros(dim, dim);
            }

            for (int i = 0; i < pointsNum; ++i) {
                int cl = labels[i];
                if (activeClusters.find(cl) != activeClusters.end())
                    sigma[cl] += (points.col(i) - cls[cl].m) * (points.col(i) - cls[cl].m).t();
            }
            std::set<int>::iterator it = activeClusters.begin();
            while (it != activeClusters.end()) {
                sigma[*it] /= cls[*it].card;
                uvec ind(dim - 1);
                for (int i = 0; i < dim - 1; ++i) ind[i] = i + 1;
                double EMin = INFINITY;
                int bestAxis;
                vec bestCoeffs;
                double bestVar;
                try {
                    for (int i = 0; i < dim; ++i) {
                        cls[*it].sigmaL[i] = chol(sigma[*it].submat(ind, ind), "lower");
                        double det = 1.0;
                        for (int j = 0; j < dim - 1; ++j) det *= cls[*it].sigmaL[i].at(j, j);

                        cls[*it].AL[i] = chol(cls[*it].A[i], "lower");
                        vec y = solve(trimatl(cls[*it].AL[i]), cls[*it].b[i]);
                        vec x = solve(trimatu(cls[*it].AL[i].t()), y);
                        double var = cls[*it].sumOfSquares[i];
                        var -= ((vec)(2.0 * x.t() * cls[*it].b[i]))[0];
                        var += ((mat)(x.t() * cls[*it].A[i] * x)).at(0, 0);
                        var /= cls[*it].card;
                        double E = cls[*it].prob * (-log(cls[*it].prob) + (0.5 * ((dim * log(2.0 * M_PI * M_E)) + log(det * det * var))));
                        if (!is_finite(E)) throw 0;
                        if (E < EMin) {
                            EMin = E; bestAxis = i; bestCoeffs = x; bestVar = var;
                        }
                    }
                    cls[*it].axis = bestAxis;
                    cls[*it].coeffs = bestCoeffs;
                    cls[*it].var = bestVar;
                    cls[*it].E = EMin;
                    ++it;
                } catch (...) {
                    std::set<int>::iterator tmp = it;
                    ++tmp;
                    inactiveClusters.insert(*it);
                    activeClusters.erase(it);
                    it = tmp;
                }
            }
        }

        void RemoveCluster(int cluster) {
            for (int i = 0; i < pointsNum; ++i) {
                if (labels[i] == cluster) {
                    double EMin = INFINITY;
                    int bestCl;
                    SClusterData bestData(*this);
                    SClusterData data(*this);
                    for (std::set<int>::iterator it = activeClusters.begin(); it != activeClusters.end(); ++it) {
                        if (*it != cluster) {
                            data = cls[*it];
                            data.Update(i);
                            if (data.E < EMin) {
                                EMin = data.E;
                                bestCl = *it;
                                bestData = data;
                            }
                        }
                    }
                    labels[i] = bestCl;
                    cls[cluster] = bestData;
                }
            }
        }

        void UpdateResult() {
            Rcpp::List cardL(maxClusters);
            Rcpp::List probL(maxClusters);
            Rcpp::List mL(maxClusters);
            Rcpp::List covL(maxClusters);
            Rcpp::List coeffsL(maxClusters);
            Rcpp::List dirsL(maxClusters);
            Rcpp::List EL(maxClusters);
            for (int i = 0; i < maxClusters; ++i) {
                if (activeClusters.find(i) != activeClusters.end()) {
                    uvec ind(dim - 1);
                    for (int j = 0; j < cls[i].axis; ++j) ind(j) = j;
                    for (int j = cls[i].axis + 1; j < dim; ++j) ind(j - 1) = j;

                    sigma[i](cls[i].axis, span::all).fill(0.0);
                    sigma[i](span::all, cls[i].axis).fill(0.0);
                    sigma[i](cls[i].axis, cls[i].axis) = cls[i].var;
                    sigma[i].submat(ind, ind) = cls[i].sigmaL[cls[i].axis] * cls[i].sigmaL[cls[i].axis].t();

                    cardL[i] = cls[i].card;
                    probL[i] = cls[i].prob;
                    mL[i] = cls[i].m;
                    covL[i] = sigma[i];
                    coeffsL[i] = cls[i].coeffs;
                    dirsL[i] = cls[i].axis;
                    EL[i] = cls[i].E;
                } else {
                    cardL[i] = R_NilValue;
                    probL[i] = R_NilValue;
                    mL[i] = R_NilValue;
                    covL[i] = R_NilValue;
                    coeffsL[i] = R_NilValue;
                    dirsL[i] = R_NilValue;
                    EL[i] = R_NilValue;
                }
            }
            if (interactive)
                res.push_back(Rcpp::List::create(
                    Rcpp::Named("number_of_clusters", (int)activeClusters.size()),
                    Rcpp::Named("labels", labels),
                    Rcpp::Named("cardinalities", cardL),
                    Rcpp::Named("probabilities", probL),
                    Rcpp::Named("means", mL),
                    Rcpp::Named("covariances", covL),
                    Rcpp::Named("coefficients", coeffsL),
                    Rcpp::Named("directions", dirsL),
                    Rcpp::Named("cost", EL),
                    Rcpp::Named("cost_total", ENew)
                ));
            else
                res = Rcpp::List::create(
                    Rcpp::Named("number_of_clusters", (int)activeClusters.size()),
                    Rcpp::Named("labels", labels),
                    Rcpp::Named("cardinalities", cardL),
                    Rcpp::Named("probabilities", probL),
                    Rcpp::Named("means", mL),
                    Rcpp::Named("covariances", covL),
                    Rcpp::Named("coefficients", coeffsL),
                    Rcpp::Named("directions", dirsL),
                    Rcpp::Named("cost", EL),
                    Rcpp::Named("cost_total", ENew),
                    Rcpp::Named("number_of_iterations", numberOfIterations)
                );
        }
};

//   -*-   -*-   -*-

#endif
