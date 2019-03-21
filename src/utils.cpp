#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

//' @useDynLib combinIT, .registration = TRUE
//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
arma::mat res(arma::mat x) {
  int nr = as_scalar(x.n_rows);
  int nc = as_scalar(x.n_cols);
  arma::vec RowMean = arma::mean(x,1);
  arma::vec ColMean = arma::trans(arma::mean(x,0));
  double Mean = as_scalar(arma::accu(x)/(nc*nr));
  arma::mat y(nr,nc);
  for(int i=0; i<nr;i++)
  {
    for(int j=0;j< nc;j++)
    {
      y(i,j) = x(i,j)-RowMean(i)-ColMean(j)+Mean; 
    }
  }
  return y;
}

//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
  float Bfc(arma::mat x,int bl, int tr,int p) {
  mat y(bl,tr),yt1(bl,tr),yt2(bl,tr); 
  y = res(x);
  yt1 = y.t() * y;
  yt2 = yt1*yt1;
  float Boik = arma::trace(yt1)*arma::trace(yt1) / (p*trace(yt2));
  return Boik;
}

//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
arma::vec Bfsim(int nsim,int bl, int tr,int p){
  mat sam(bl,tr);
  vec out(nsim);
  for(int i=0;i<nsim;i++)
  {
    sam.randn(bl,tr);
    out(i)=Bfc(sam,bl,tr,p);
  }
  return out;
}
//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
double picf(arma::vec y,arma::mat kp,float c0){
  vec z= kp * y;
  for(unsigned int i=0;i<kp.n_rows;i++)
    z(i)=fabs(z(i));
    vec s0=median(z,0);
    uvec ids = find(z <= (5*s0(0)) );
    vec PSE=median(z.elem(ids),0);
    vec PIC=max(z,0)/PSE(0);
    return PIC(0);
}
//' Module  function
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
arma::vec PICfsim(int nsim,arma::mat kp, float c0, int n){
  vec sam(n);
  vec out(nsim);
  for(int i=0;i<nsim;i++)
  {
    sam.randn(n);
    out(i)=picf(sam,kp,c0);
  }
  return out;
}

//' @importFrom Rcpp sourceCpp
//'
// [[Rcpp::export]]
double C0(arma::mat kp, int n,int nc0){
  vec sim(nc0);
  vec norsam(n);
  for(int i=0;i<nc0;i++)
  {
    vec temp=kp*norsam.randn(n);
    for(unsigned int j=0;j<temp.n_rows;j++)
    temp(j) = fabs(temp(j));
    vec me=median(temp,0);
    sim(i)=me(0);
  }
  vec out = mean(sim,0);
  return out(0);
}

//' @importFrom Rcpp sourceCpp
//' @useDynLib combinIT
//' 
// [[Rcpp::export]]
double piephoC(arma::mat x,int bl, int tr) {
  mat Res(bl,tr),Res2(bl,tr);
  Res =res(x);
  Res2 =Res%Res;
  vec RowSum = arma::sum(Res2,1);
  vec delta = (bl*(bl-1)*RowSum-sum(RowSum));
  double h1 = 0;
  for(int i=0;i< (bl-1);i++)
    for(int j=i+1; j < bl ; j++)
      h1 += as_scalar(delta(i)*delta(j));
    double U = as_scalar(2*bl*h1/((bl-1)*pow(sum(delta),2)));
    double piepho = -(tr-1)*(bl-1)*(bl-2)*log(U)/2;
    return piepho;
}

//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
arma::vec Piephosim(int nsim,int bl, int tr){
  mat sam(bl,tr);
  vec out(nsim);
  for(int i=0;i<nsim;i++)
  {
    sam.randn(bl,tr);
    out(i)=piephoC(sam,bl,tr);
  }
  return out;
}

//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
arma::umat mycombn(int n, int k) {
  double n_subsets = Rf_choose(n, k);
  umat out = zeros<umat>(k, n_subsets);
  uvec a = linspace<uvec>(1, k, k);  
  out.col(0) = a;
  int m = 0;
  int h = k;
  uvec j;
  for(long long i = 1; i < n_subsets; i++){
    if(m < (n - h)){  
      h = 1;
      m = a(k - 1);
      j = linspace<uvec>(1, 1, 1);
    }
    else{
      m = a(k - h - 1);
      ++h;
      j = linspace<uvec>(1, h, h);
    }
    a.elem(k - h - 1 + j) = m + j;
    out.col(i) = a;
  }
  return(out);
}


//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
Rcpp::LogicalVector logical_index(Rcpp::IntegerVector idx, R_xlen_t n) {
  bool invert = false; 
  Rcpp::LogicalVector result(n, false);
  for (R_xlen_t i = 0; i < idx.size(); i++) {
    if (!invert && idx[i] < 0) invert = true;
    result[std::abs(idx[i])] = true;
  }
  
  if (!invert) return result;
  return !result;
}
//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
NumericVector  Subset(Rcpp::NumericVector x, Rcpp::IntegerVector idx) {
    return x[logical_index(idx, x.size())];
  }

//' @importFrom Rcpp sourceCpp
//' @useDynLib combinIT
//' 
// [[Rcpp::export]]
double kkf_C(arma::mat x,int bl, int tr)
{
   mat yb1,yb2,yb11,yb22;
   int count = 0;
   NumericVector Sel(bl);
   NumericVector rows_nsel; 
   double fvalues=0;
   vec pvalues(1);
   double rss1,rss2;
   for(int i=0;i<bl;i++)
     Sel[i]=i;
   for(int i=1; i< floor(bl/2);i++){
    umat ind = mycombn(bl,i+1) -1 ;
    int Nsplit = ind.n_cols;
    if((bl/2)== (i)) Nsplit = Nsplit/2;
     for(int j=0;j<Nsplit;j++){
        yb1 = x.rows(ind.col(j));
        IntegerVector nSel = as<IntegerVector>(wrap(ind.col(j)));
        rows_nsel =  Subset(Sel,-1*nSel);
        yb2 = x.rows(as<arma::uvec>(rows_nsel)); 
        yb1 = res(yb1);
        yb11 = yb1%yb1;
        rss1 = accu(yb11);
        yb2 = res(yb2);
        yb22 = yb2%yb2;
        rss2 = accu(yb22);
        double dfn = as_scalar((tr-1)*(i+1-1));
        double dfd = as_scalar((bl-(i+1)-1)*(tr-1));
        fvalues = (rss1*(bl-(i+1)-1))/(rss2*(i+1-1));
        if(fvalues<1) fvalues = 1/fvalues;
        pvalues[count]= 1-R::pf(fvalues,dfn,dfd,1,0)+R::pf(1/fvalues,dfn,dfd,1,0);
        count ++; 
        pvalues.resize(pvalues.size()+1);
     }
   }
   double KKSA = min(pvalues.elem(find(pvalues>1e-20)));
   return KKSA;
 }

//' @importFrom Rcpp sourceCpp
//' @useDynLib combinIT
//' 
// [[Rcpp::export]]
NumericVector KKsim(int nsim,int bl, int tr){
  mat sam(bl,tr);
  NumericVector out(nsim);
  for(int i=0;i< nsim;i++)
  {
    sam.randn(bl,tr);
    out(i) = kkf_C(sam,bl,tr);
  }
  return out;
}

/*
//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
double hf_C(arma::mat x,int bl, int tr)
{
  double sse =as_scalar(sum(pow(res(x),2))); 
  vec hvalues();
  mat yb1,yb2;
  int count = 0;
  for(int i=1; i< floor(bl/2);i++){
    umat ind = mycombn(bl,i+1) -1 ;
    int Nsplit = ind.n_cols;
    if((bl/2)== (i)) Nsplit = Nsplit/2;
    for(int j=0;j<Nsplit;j++){
      count ++;
      yb1 = x.rows(ind.col(j));               // check
      yb2 = x.rows(ind.col(j));
      double rss1 = as_scalar(sum(pow(res(yb1),2)));
      double rss2 = as_scalar(sum(pow(res(yb2),2)));
      double sse7  = rss1+rss2;
      hvalues(count) = (sse-sse7)*(bl-2)/sse7;
      
    }
  }
  for(int i=0; i< bl;i++){
    count ++;
    yb1 = x.rows(-i);
    double rss7 = as_scalar(sum(pow(res(yb1),2)));
    hvalues(count) = (sse-sse7)*(bl-2)/sse7;
  }
  
  double fmax = max(hvalues,0);
  return(fmax);
  
}

//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
arma::vec hfsim(int nsim,int bl, int tr){
  mat sam(bl,tr);
  vec out(nsim);
  for(int i=0;i<nsim;i++)
  {
    sam.randn(bl,tr);
    out(i)=hf_C(sam,bl,tr);
  }
  return out;
}
*/ 
