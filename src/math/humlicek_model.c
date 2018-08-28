/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void humlicek_model(double x, double m, double S, double A, double eta,
		    double *f,
		    double *dfdm, double *dfds, double *dfdA, double *dfdeta) {
  double X = (x - m)/(sqrt(2.0) * S);
  double Y = fabs(eta * S * sqrt(2.0));

  double t[6] = {.314240376,.947788391,1.59768264,2.27950708,3.02063703,3.8897249};
  double c[6] = {1.01172805,-.75197147,1.2557727e-2,1.00220082e-2,-2.42068135e-4,
		 5.00848061e-7};
  double s[6] = {1.393237,.231152406,-.155351466,6.21836624e-3,9.19082986e-5,
		 -6.27525958e-7};
  double Vreal = 0.0;
  double Vimag = 0.0;
  double wr = 0.0;
  double wi = 0.0;
  double y1 = Y + 1.5;
  double y2 = y1 * y1;
  int i;
  double y3,r,r2,d,d1,d2,d3,d4;
  double dVdX,dVdY;
  if ((Y < 0.85)&&(fabs(X) > 18.1 * Y + 1.65)) {
    if (fabs(X) < 12.0) {
      wr = exp(-1.0 * X * X);
    }
    y3 = Y + 3.0;
    for (i = 0; i < 6; i++) {
      r = X - t[i];
      r2= r * r;
      d = 1.0 / (r2 + y2);
      d1= y1 * d;
      d2= r * d;
      wr += Y * (c[i] * (r * d2 - 1.5 * d1) + s[i] * y3 * d2) / (r2 + 2.25);
      r  = X + t[i];
      r2 = r * r;
      d  = 1.0 / (r2 + y2);
      d3 = y1 * d;
      d4 = r * d;
      wr += Y * (c[i] * (r * d4 - 1.5 * d3) - s[i] * y3 * d4) / (r2 + 2.25);
      wi += c[i] * (d2 + d4) + s[i] * (d1 - d3);
    }
  }
  else {
    for (i = 0; i < 6; i++) {
      r = X - t[i];
      r2= r * r;
      d = 1.0 / (r2 + y2);
      d1= y1 * d;
      d2= r * d;
      r  = X + t[i];
      r2 = r * r;
      d  = 1.0 / (r2 + y2);
      d3 = y1 * d;
      d4 = r * d;
      wr += c[i] * (d1 + d3) - s[i] * (d2 - d4);
      wi += c[i] * (d2 + d4) + s[i] * (d1 - d3);
    }
  }

  Vreal = wr;
  Vimag = wi;
  dVdX = 2.0 * (Y * Vimag - X * Vreal);
  dVdY = 2.0 * (X * Vimag + Y * Vreal) - 2.0 * M_1_PI;
  
  *dfdA   = Vreal / sqrt(2.0 * M_PI * ((1 + eta) * S));
  *f      = A * *dfdA * sqrt(2.0 * M_PI * ((1 + eta) * S)) ;
  *dfdeta = -dVdY;
  *dfdm   = dVdX / (sqrt(2.0) * S);
  *dfds   = -dVdX * (m - x)/(sqrt(2) * pow((1 + eta) * S,2)) +
    eta * *dfdeta / ( S);

}

double humlicek_line(double x, lines *L, int i) {
  double f,dfdA,dfdm,dfds,dfdeta;
  humlicek_model(x,L->m[i],L->s[i],L->F[i]  * pow(L->s[i] * sqrt(2.0 * M_PI),-1)
		 ,L->eta[i],
		 &f,&dfdm,&dfds,&dfdA,&dfdeta);
  return(f);
}
