// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// BayesMVP_internal
int BayesMVP_internal(const std::string& dataFile, const std::string& mrfGFile, const std::string& blockFile, const std::string& structureGraphFile, const std::string& hyperParFile, const std::string& outFilePath, unsigned int nIter, unsigned int burnin, unsigned int nChains, const std::string& covariancePrior, const std::string& gammaPrior, const std::string& gammaSampler, const std::string& gammaInit, const std::string& wSampler, const std::string& betaPrior, const int maxThreads, const int tick, bool output_gamma, bool output_beta, bool output_Gy, bool output_sigmaRho, bool output_pi, bool output_tail, bool output_model_size, bool output_CPO, bool output_model_visit);
RcppExport SEXP _BayesMVP_BayesMVP_internal(SEXP dataFileSEXP, SEXP mrfGFileSEXP, SEXP blockFileSEXP, SEXP structureGraphFileSEXP, SEXP hyperParFileSEXP, SEXP outFilePathSEXP, SEXP nIterSEXP, SEXP burninSEXP, SEXP nChainsSEXP, SEXP covariancePriorSEXP, SEXP gammaPriorSEXP, SEXP gammaSamplerSEXP, SEXP gammaInitSEXP, SEXP wSamplerSEXP, SEXP betaPriorSEXP, SEXP maxThreadsSEXP, SEXP tickSEXP, SEXP output_gammaSEXP, SEXP output_betaSEXP, SEXP output_GySEXP, SEXP output_sigmaRhoSEXP, SEXP output_piSEXP, SEXP output_tailSEXP, SEXP output_model_sizeSEXP, SEXP output_CPOSEXP, SEXP output_model_visitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type dataFile(dataFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type mrfGFile(mrfGFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type blockFile(blockFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type structureGraphFile(structureGraphFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type hyperParFile(hyperParFileSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type outFilePath(outFilePathSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nIter(nIterSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nChains(nChainsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type covariancePrior(covariancePriorSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type gammaPrior(gammaPriorSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type gammaSampler(gammaSamplerSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type gammaInit(gammaInitSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type wSampler(wSamplerSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type betaPrior(betaPriorSEXP);
    Rcpp::traits::input_parameter< const int >::type maxThreads(maxThreadsSEXP);
    Rcpp::traits::input_parameter< const int >::type tick(tickSEXP);
    Rcpp::traits::input_parameter< bool >::type output_gamma(output_gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type output_beta(output_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type output_Gy(output_GySEXP);
    Rcpp::traits::input_parameter< bool >::type output_sigmaRho(output_sigmaRhoSEXP);
    Rcpp::traits::input_parameter< bool >::type output_pi(output_piSEXP);
    Rcpp::traits::input_parameter< bool >::type output_tail(output_tailSEXP);
    Rcpp::traits::input_parameter< bool >::type output_model_size(output_model_sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type output_CPO(output_CPOSEXP);
    Rcpp::traits::input_parameter< bool >::type output_model_visit(output_model_visitSEXP);
    rcpp_result_gen = Rcpp::wrap(BayesMVP_internal(dataFile, mrfGFile, blockFile, structureGraphFile, hyperParFile, outFilePath, nIter, burnin, nChains, covariancePrior, gammaPrior, gammaSampler, gammaInit, wSampler, betaPrior, maxThreads, tick, output_gamma, output_beta, output_Gy, output_sigmaRho, output_pi, output_tail, output_model_size, output_CPO, output_model_visit));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesMVP_BayesMVP_internal", (DL_FUNC) &_BayesMVP_BayesMVP_internal, 26},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesMVP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
