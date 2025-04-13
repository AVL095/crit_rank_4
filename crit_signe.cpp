#include <cmath>
#include <fstream>

using namespace std;

// Got this from http://www.math.ucla.edu/~tom/distributions/binomial.html
double betinc(double x,double a,double b,double eps) {
  double a0,b0,a1,b1,m9,a2,c9;
  a0 = 0; b0 = 1; a1 = 1; b1 = 1; m9 = 0; a2 = 0;
  while (abs((a1 - a2) / a1) > eps) {
    a2 = a1;
    c9 = -(a + m9) * (a + b + m9) * x / (a + 2. * m9) / (a + 2. * m9 + 1.);
    a0 = a1 + c9 * a0;
    b0 = b1 + c9 * b0;
    m9 = m9 + 1.;
    c9 = m9 * (b - m9) * x / (a + 2. * m9 - 1) / (a + 2. * m9);
    a1 = a0 + c9 * a1;
    b1 = b0 + c9 * b1;
    a0 = a0 / b1;
    b0 = b0 / b1;
    a1 = a1 / b1;
    b1 = 1.;
  }
  return a1 / a;
}
///////////////////
double binomial_cdf(double x, int n, double p) {
    double betacdf;
    double eps = 1e-10;

    if (x < 0)    return 0.0;
    if (x >= n)    return 1.0;
    if (p < 0.0 || p > 1.0 || n <= 0)  return numeric_limits<double>::quiet_NaN();

    x = floor(x);
    double z = p;
    double a = x + 1.;
    double b = n - x;
    double s = a + b;
    double bt = exp(lgamma(s) - lgamma(b) - lgamma(a) + a * log(z) + b * log(1 - z));
    if (z < (a + 1.) / (s + 2.))
        betacdf = bt * betinc(z, a, b, eps);
    else
        betacdf = 1.0 - bt * betinc(1.0 - z, b, a, eps);
    return round((1.0 - betacdf) * (1.0 / eps)) / (1.0 / eps);
}

////////////////////////////////////////////////////////////////////////////////

int main() {
   int n,i;
   double *pw,*w,s;
    n=45;

   w=new double[n+1];
   pw=new double[n+1];
   for(i=0;i<=n;i++) {
      w[i]=double(i);
      pw[i]=binomial_cdf(w[i],n,0.5);
   }
   s=0;
    for(i=0;i<=n;i++) {
     s+=pw[i]/pow(2,n);
    } 

    ofstream out("signe.out");
    for(i=0;i<=n;i++) out<<w[i]<<";"<<pw[i]<<endl;

   out.close();
   delete []w,pw;
   return 0;
}

