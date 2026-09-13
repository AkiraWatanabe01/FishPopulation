#ifndef _FISHCODE_H_
#define _FISHCODE_H_
 
double f_Nl(double x1,  double m, double f);
double f_Nh(double x1, double x2, double m, double f1, double f2);
double f_B(double x, double ma, double w);
double f_W(double x, double w);
double f_Cn(double x, double mo, double f);
double f_Cw(double x, double w, double mo, double f);
double Ricker(double ssb, double alp, double bet);
double Hockey(double ssb, double alp, double Rinf);
double Beverton(double ssb, double alp, double bet);
double f_NCw(double c, double w, double mo, double f);
double f_NCn(double c, double mo, double f);
double maxi(double *comp, int clen);
void mult_TAC(double tac, double *rnd, double *x, double *w, double *f, int age, double mor);
double mult_ABC(double ssb, double ssb_01, double ssb_msy, double ssb_06);
void CalculateAbundanceNextyear(double abundance_next_year, double abundance, double natural_mortality, double fishing_mortality);
double CalculateRPS(int type, double ssb, int fyr);
double CalculateCoefficientForABC(double ssb, double ssb_01, double ssb_msy, double ssb_06);

#endif // _FISHCODE_H_
