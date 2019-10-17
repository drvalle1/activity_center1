#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function helps with multinomial draws
int cat1(double value, NumericVector prob) {
  int res=prob.length();
  double probcum = 0;
  
  for (int i = 0; i < prob.length(); i++) {
    probcum = probcum + prob(i);
    if (value < probcum) {
      res = i;
      break;
    }
  }
  return res;
}

//' This function samples z's from a categorical distribution
// [[Rcpp::export]]
IntegerVector rmultinom1(NumericMatrix prob, NumericVector randu) {
  
  IntegerVector z(prob.nrow());

  for(int i=0; i<prob.nrow();i++){
    z[i]=cat1(randu[i],prob(i,_));
  }
  return z;
}

//' This function calculates distance between two sets of coordinates
// [[Rcpp::export]]
NumericMatrix GetDistance(NumericMatrix AcCoord,NumericMatrix GridCoord, int Ngrid, int Nac) {
  
  NumericMatrix res(Ngrid,Nac);
  double x2;
  double y2;
  for(int i=0; i<Ngrid;i++){
    for(int j=0; j<Nac; j++){
      x2=pow(GridCoord(i,0)-AcCoord(j,0),2.0);
      y2=pow(GridCoord(i,1)-AcCoord(j,1),2.0);
      res(i,j)=sqrt(x2+y2);
    }
  }
  return res;
}
