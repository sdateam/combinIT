// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Seq
IntegerVector Seq(int a, int b);
RcppExport SEXP _combinIT_Seq(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(Seq(a, b));
    return rcpp_result_gen;
END_RCPP
}
// mod
NumericVector mod(int K, int W);
RcppExport SEXP _combinIT_mod(SEXP KSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(mod(K, W));
    return rcpp_result_gen;
END_RCPP
}
// Csij
double Csij(int i, int j, int K, int L, NumericMatrix B);
RcppExport SEXP _combinIT_Csij(SEXP iSEXP, SEXP jSEXP, SEXP KSEXP, SEXP LSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(Csij(i, j, K, L, B));
    return rcpp_result_gen;
END_RCPP
}
// SS
NumericMatrix SS(int K, int L, NumericMatrix B, int d);
RcppExport SEXP _combinIT_SS(SEXP KSEXP, SEXP LSEXP, SEXP BSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(SS(K, L, B, d));
    return rcpp_result_gen;
END_RCPP
}
// Cofmat
NumericMatrix Cofmat(int d, int L, NumericVector cx);
RcppExport SEXP _combinIT_Cofmat(SEXP dSEXP, SEXP LSEXP, SEXP cxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cx(cxSEXP);
    rcpp_result_gen = Rcpp::wrap(Cofmat(d, L, cx));
    return rcpp_result_gen;
END_RCPP
}
// Gram
NumericMatrix Gram(int K, int L, NumericMatrix A, int d);
RcppExport SEXP _combinIT_Gram(SEXP KSEXP, SEXP LSEXP, SEXP ASEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type L(LSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(Gram(K, L, A, d));
    return rcpp_result_gen;
END_RCPP
}
// res
arma::mat res(arma::mat x);
RcppExport SEXP _combinIT_res(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(res(x));
    return rcpp_result_gen;
END_RCPP
}
// Bfc
float Bfc(arma::mat x, int bl, int tr, int p);
RcppExport SEXP _combinIT_Bfc(SEXP xSEXP, SEXP blSEXP, SEXP trSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type bl(blSEXP);
    Rcpp::traits::input_parameter< int >::type tr(trSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(Bfc(x, bl, tr, p));
    return rcpp_result_gen;
END_RCPP
}
// Bfsim
arma::vec Bfsim(int nsim, int bl, int tr, int p);
RcppExport SEXP _combinIT_Bfsim(SEXP nsimSEXP, SEXP blSEXP, SEXP trSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type bl(blSEXP);
    Rcpp::traits::input_parameter< int >::type tr(trSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(Bfsim(nsim, bl, tr, p));
    return rcpp_result_gen;
END_RCPP
}
// picf
double picf(arma::vec y, arma::mat kp, float c0);
RcppExport SEXP _combinIT_picf(SEXP ySEXP, SEXP kpSEXP, SEXP c0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type kp(kpSEXP);
    Rcpp::traits::input_parameter< float >::type c0(c0SEXP);
    rcpp_result_gen = Rcpp::wrap(picf(y, kp, c0));
    return rcpp_result_gen;
END_RCPP
}
// PICfsim
arma::vec PICfsim(int nsim, arma::mat kp, float c0, int n);
RcppExport SEXP _combinIT_PICfsim(SEXP nsimSEXP, SEXP kpSEXP, SEXP c0SEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type kp(kpSEXP);
    Rcpp::traits::input_parameter< float >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(PICfsim(nsim, kp, c0, n));
    return rcpp_result_gen;
END_RCPP
}
// C0
double C0(arma::mat kp, int n, int nc0);
RcppExport SEXP _combinIT_C0(SEXP kpSEXP, SEXP nSEXP, SEXP nc0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type kp(kpSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type nc0(nc0SEXP);
    rcpp_result_gen = Rcpp::wrap(C0(kp, n, nc0));
    return rcpp_result_gen;
END_RCPP
}
// piephoC
double piephoC(arma::mat x, int bl, int tr);
RcppExport SEXP _combinIT_piephoC(SEXP xSEXP, SEXP blSEXP, SEXP trSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type bl(blSEXP);
    Rcpp::traits::input_parameter< int >::type tr(trSEXP);
    rcpp_result_gen = Rcpp::wrap(piephoC(x, bl, tr));
    return rcpp_result_gen;
END_RCPP
}
// Piephosim
arma::vec Piephosim(int nsim, int bl, int tr);
RcppExport SEXP _combinIT_Piephosim(SEXP nsimSEXP, SEXP blSEXP, SEXP trSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type bl(blSEXP);
    Rcpp::traits::input_parameter< int >::type tr(trSEXP);
    rcpp_result_gen = Rcpp::wrap(Piephosim(nsim, bl, tr));
    return rcpp_result_gen;
END_RCPP
}
// mycombn
arma::umat mycombn(int n, int k);
RcppExport SEXP _combinIT_mycombn(SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(mycombn(n, k));
    return rcpp_result_gen;
END_RCPP
}
// logical_index
Rcpp::LogicalVector logical_index(Rcpp::IntegerVector idx, R_xlen_t n);
RcppExport SEXP _combinIT_logical_index(SEXP idxSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< R_xlen_t >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(logical_index(idx, n));
    return rcpp_result_gen;
END_RCPP
}
// Subset
NumericVector Subset(Rcpp::NumericVector x, Rcpp::IntegerVector idx);
RcppExport SEXP _combinIT_Subset(SEXP xSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(Subset(x, idx));
    return rcpp_result_gen;
END_RCPP
}
// kkf_C
double kkf_C(arma::mat x, int bl, int tr);
RcppExport SEXP _combinIT_kkf_C(SEXP xSEXP, SEXP blSEXP, SEXP trSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type bl(blSEXP);
    Rcpp::traits::input_parameter< int >::type tr(trSEXP);
    rcpp_result_gen = Rcpp::wrap(kkf_C(x, bl, tr));
    return rcpp_result_gen;
END_RCPP
}
// KKsim
NumericVector KKsim(int nsim, int bl, int tr);
RcppExport SEXP _combinIT_KKsim(SEXP nsimSEXP, SEXP blSEXP, SEXP trSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int >::type bl(blSEXP);
    Rcpp::traits::input_parameter< int >::type tr(trSEXP);
    rcpp_result_gen = Rcpp::wrap(KKsim(nsim, bl, tr));
    return rcpp_result_gen;
END_RCPP
}
// hf_C
List hf_C(arma::mat x, int bl, int tr);
RcppExport SEXP _combinIT_hf_C(SEXP xSEXP, SEXP blSEXP, SEXP trSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type bl(blSEXP);
    Rcpp::traits::input_parameter< int >::type tr(trSEXP);
    rcpp_result_gen = Rcpp::wrap(hf_C(x, bl, tr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_combinIT_Seq", (DL_FUNC) &_combinIT_Seq, 2},
    {"_combinIT_mod", (DL_FUNC) &_combinIT_mod, 2},
    {"_combinIT_Csij", (DL_FUNC) &_combinIT_Csij, 5},
    {"_combinIT_SS", (DL_FUNC) &_combinIT_SS, 4},
    {"_combinIT_Cofmat", (DL_FUNC) &_combinIT_Cofmat, 3},
    {"_combinIT_Gram", (DL_FUNC) &_combinIT_Gram, 4},
    {"_combinIT_res", (DL_FUNC) &_combinIT_res, 1},
    {"_combinIT_Bfc", (DL_FUNC) &_combinIT_Bfc, 4},
    {"_combinIT_Bfsim", (DL_FUNC) &_combinIT_Bfsim, 4},
    {"_combinIT_picf", (DL_FUNC) &_combinIT_picf, 3},
    {"_combinIT_PICfsim", (DL_FUNC) &_combinIT_PICfsim, 4},
    {"_combinIT_C0", (DL_FUNC) &_combinIT_C0, 3},
    {"_combinIT_piephoC", (DL_FUNC) &_combinIT_piephoC, 3},
    {"_combinIT_Piephosim", (DL_FUNC) &_combinIT_Piephosim, 3},
    {"_combinIT_mycombn", (DL_FUNC) &_combinIT_mycombn, 2},
    {"_combinIT_logical_index", (DL_FUNC) &_combinIT_logical_index, 2},
    {"_combinIT_Subset", (DL_FUNC) &_combinIT_Subset, 2},
    {"_combinIT_kkf_C", (DL_FUNC) &_combinIT_kkf_C, 3},
    {"_combinIT_KKsim", (DL_FUNC) &_combinIT_KKsim, 3},
    {"_combinIT_hf_C", (DL_FUNC) &_combinIT_hf_C, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_combinIT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
