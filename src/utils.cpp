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
arma::mat kkf_C(arma::mat x,int bl, int tr)
{
  //IntegerVector Nrow =  seq(2, floor(bl/2));
   //return as<NumericVector>(Nrow); 
   vec fvalues ();
   vec pvalues();
   mat yb1,yb2;
   int count = 0;
   for(int i=1; i< floor(bl/2);i++){
     umat ind = mycombn(bl,i+1) -1 ;
     int Nsplit = ind.n_cols;
    if((bl/2)== (i)) Nsplit = Nsplit/2;
     for(int j=0;j<Nsplit;j++){
       count ++;
        yb1 = x.rows(ind.col(j));
        vec colj = ind.col(j);
        IntegerVector col_j = seq(1,ind.n_cols) 
        yb2 = x.rows();
       //rss1<-sum(( t(yb1 - apply(yb1, 1, mean) + mean(yb1)) - apply(yb1, 2, mean))^2)
       //rss2<-sum(( t(yb2 - apply(yb2, 1, mean) + mean(yb2)) - apply(yb2, 2, mean))^2)
       //dfn<-(tr-1)*(i-1)
       //dfd<-(bl-i-1)*(tr-1)
       //fvalues[count]<-(rss1*(bl-i-1))/(rss2*(i-1))
       //if(fvalues[count]<1)fvalues[count]<-1/fvalues[count]
       //pvalues[count]<-1-pf(fvalues[count],dfn,dfd)+pf(1/fvalues[count],dfn,dfd)
     }
   }
   //KKSA<-min(pvalues)
     return yb2;
   
}