// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;
using namespace Rcpp;

//' calculates factor-stratified means for each column
//'
//' @param sY sparse matrix Gene by cell matrix of counts
//' @param rowSel factor The donor that each cell is from
//' @param numCells numeric Total number of cells belonging to each donor
//' @export
// [[Rcpp::export]]
NumericVector get_vars(SEXP sY,  IntegerVector rowSel, IntegerVector numCells) {

  // need to do this as SEXP, modify the slots on the fly
  S4 mat(sY);
  const arma::uvec i(( unsigned int *)INTEGER(mat.slot("i")),LENGTH(mat.slot("i")),false,true);
  const arma::ivec dims(INTEGER(mat.slot("Dim")),LENGTH(mat.slot("Dim")),false,true);
  const arma::ivec p(INTEGER(mat.slot("p")),LENGTH(mat.slot("p")),false,true);
  arma::vec Y(REAL(mat.slot("x")),LENGTH(mat.slot("x")),false,true);

  List dimnames(mat.slot("Dimnames"));
  CharacterVector geneNames(dimnames[1]);

  CharacterVector factorLevels=rowSel.attr("levels");
  int nlevels=factorLevels.size();
  CharacterVector expandedFactorLevels(nlevels+1);
  expandedFactorLevels[0]="NA";
  for(int i=0;i<nlevels;i++) { expandedFactorLevels[i+1]=factorLevels[i]; }
  const arma::ivec rs=arma::ivec(INTEGER(rowSel),LENGTH(rowSel),false,true);

  int ncols=p.size()-1;

  if(nlevels==0) { stop("colSumByFac(): supplied factor doesn't have any levels!"); }
  //arma::mat sumM(nlevels+1,ncols,arma::fill::zeros);
  NumericMatrix sumM(nlevels+1,ncols);

  // for each gene
  //#pragma omp parallel for shared(sumM)
  for(int g=0;g<ncols;g++) {
    int p0=p[g]; int p1=p[g+1];
    if(p1-p0 <1) { continue; }
    for(int j=p0;j<p1;j++) {
      int row=i[j];
      int f=rs[row];
      if(f==NA_INTEGER) {
        sumM(0,g)+=Y[j];
      } else if(f>0) {
        sumM(f,g)+=Y[j]/numCells[f-1];
      }
    }
  }


  colnames(sumM) = geneNames;
  rownames(sumM) = expandedFactorLevels;

  // calculate the variance across donors (rows)
  int nrow = sumM.nrow(), ncol = sumM.ncol();
  NumericVector out(ncol);

  for (int j = 0; j < ncol; j++) {
    double mean = 0;
    double M2 = 0;
    int n;
    double delta, xx;
    for (int i = 1; i < nrow; i++) {
      n = i;
      xx = sumM(i,j);
      delta = xx - mean;
      mean += delta/n;
      M2 = M2 + delta*(xx-mean);
    }
    out(j) = M2/(n-1);
  }
  return out;

}






