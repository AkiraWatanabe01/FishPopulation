// FishCode.c
// Cohort Analysis
#include<math.h>
#include<stdio.h>

// ======= oa ========
// == Ricker (1976 to 1987) ==
#define RI_ALPHA_HIGH (double)(56.67e0)    // (kg)^{-1}
#define RI_BETA_HIGH (double)(1e-4)       // (kilo ton)^{-1}
// == RI (1988 to 2019) ==
#define RI_ALPHA_LOW (double)(44.71e0)     // (kg)^{-1}
#define RI_BETA_LOW (double)(4.939e-4)    // (kilo ton)^{-1}
// ====== pia =======
// == BH (1976 to 1987) ==
#define BH_BETA_HIGH (double)(2.274e-4)       // (kilo ton)^{-1}
// == Bh (1988 to 2019) ==
#define BH_BETA_LOW (double)(2.563e-3)    // (kilo ton)^{-1}
// ======= ob ========
// == Ricker (1976 to 1987) ==
#define RI_ALPHA_HIGH_B (double)(56.67e0)    // (kg)^{-1}
#define RI_BETA_HIGH_B (double)(1e-4)       // (kilo ton)^{-1}
// == Hockey Stick (1988 to 1991) ==
#define HO_ALPHA_B (double)(1.86e0)     // (kg)^{-1}
#define HO_BETA_B (double)(5e4)    // (kilo ton)^{-1}
// == Ricker (1992 to 2019) ==
#define RI_ALPHA_MID_B (double)(34.3e0)    // (kg)^{-1}
#define RI_BETA_MID_B (double)(2.361e-4)       // (kilo ton)^{-1}
// ======== pib =======
// == BH (1976 to 1987) ==
#define BH_BETA_HIGH_B (double)(2.274e-4)       // (kilo ton)^{-1}
// == Bh (1988 to 1991) ==
#define BH_BETA_LOW_B (double)(4e-5)    // low regime
// == Bh (1992 to 2019) ==
#define BH_BETA_MID_B (double)(2.728e-4)    // middle regime

double f_Nl(double x1,  double m, double f) {
	return (x1 * exp(-1 * m - f));
}

double f_Nh(double x1, double x2, double m, double f1, double f2) {
	return (x1 * exp(-1 * m - f1) + x2 * exp(-1 * m - f2));
}

double f_B(double x, double ma, double w) {  // million_bi * (dimensionless) * g/bi
	return (x * w * ma * 1e-3);
}

double f_W(double x, double w) {  // million_bi * (dimensionless) * g/bi
	return (x * w * 1e-3);
}

double f_Cn(double x, double mo, double f) {   // million_bi * g/bi * () * ()
	return (x * (1 - exp(-1 * f)) * exp(-0.5 * mo));
}

double f_Cw(double x, double w, double mo, double f) {   // million_bi * g/bi * () * ()
	return (x * w * (1-exp(-1 * f)) * exp(-0.5 * mo) * 1e-3);
}

// ===== Ricker model =====
double Ricker(double ssb, double alpha, double beta) {
	return (alpha * exp(-1 * beta * ssb));
}

// ===== *Hockey Stick model =====
double Hockey(double ssb, double alpha, double Rinf) {
	return (Rinf * (1 - exp(-1 * alpha * ssb / Rinf)));  
}

// ===== Beverton Holt model =====
double Beverton(double ssb, double alpha, double beta){
	return alpha / (1 + beta * ssb);  
}

double CalculateRPS(int type, double ssb, int fyr) {
	double rps;

	switch(type){
	case 1:
        if (fyr < 1988) {
            return Ricker(ssb_true, RI_ALPHA_HIGH, RI_BETA_HIGH);
	    }
		rps = Ricker(ssb_true,RI_ALPHA_LOW,RI_BETA_LOW);
        break;
	case 2:
        if (fyr < 1988) {
        	return Beverton(ssb_true, alp[ipD+1], BH_BETA_HIGH);
	    }
        rps = Beverton(ssb_true, alp[ipD+1], BH_BETA_LOW);
	    break;
	case 3:
        if (fyr < 1988) {
            return Ricker(ssb_true, RI_ALPHA_HIGH_B, RI_BETA_HIGH_B);
	    }
		if ((fyr >= 1988) && (fyr < 1992)) {
            return Hockey(ssb_true, HO_ALPHA_B, HO_BETA_B) / ssb_true;
        }
		rps = Ricker(ssb_true, RI_ALPHA_MID_B, RI_BETA_MID_B);
	    break;
	case 4:
        if (fyr < 1988) {
            return Beverton(ssb_true, alp[ipD+1], BH_BETA_HIGH_B);
	    }
		if ((fyr >= 1988) && (fyr < 1992)) {
            return Beverton(ssb_true, alp[ipD+1], BH_BETA_LOW_B);
        }
		rps=Beverton(ssb_true, alp[ipD+1], BH_BETA_MID_B);
	    break;
	default:
		rps = 0e0;
        break;
    }

	return rps; 
}

double f_NCw(double c, double w, double mo, double f){   // million_bi * g/bi * () * ()
	return (1e+3 * c * exp(0.5 * mo) / (w * (1 -exp(-1 * f))));
}

double f_NCn(double c, double mo, double f) {   // million_bi * g/bi * () * ()
	return (c * exp(0.5 * mo) / (1 -exp(-1 * f)));
}

double maxi(double *comp, int clen) {
	double maxi = -1e0;
	int i;

	for (i = 0; i < clen; i++){
		if (maxi < comp[i]) {
			maxi=comp[i];
   		}
	}
	return maxi;
}

double CalculateCoefficientForABC(double ssb, double ssb_01, double ssb_msy, double ssb_06) {     
	if (ssb <= ssb_01){
		return 0e0;
	} 

	if ((ssb > ssb_01) && (ssb <= ssb_06)){
		return (0.8e0 * (ssb - ssb_01) / (ssb_06 - ssb_01));
	}

	return 0.8e0;
} 

void mult_TAC(double tac, double *rnd, double *x, double *w, double *f, int age, double mor){
	double nmb5=0.01, presum=0e0;
	double f_asm[age], c_asm[age];
    double c_asm_sum;
	double dif;
	short presgn, k;
	short sgn=0;   
		
	do {
		presgn=sgn;	
        c_asm_sum=0e0;

		for (k = 0; k < age; k++){
			f_asm[k] = rnd[k] * nmb5;
            if (f_asm[k] >= 3.5e0) {
				f_asm[k]=3.5e0;
			}
            
			c_asm[k] = f_Cw(x[k], w[k], mor, f_asm[k]);
            c_asm_sum += c_asm[k];
	    }

		dif=c_asm_sum-presum;
		if (dif < 0) {
			dif = -1e0 * dif;
		}
        
		if (dif<=1e-3){
			printf("break dif\n");
			break;
		}	   
        presum=c_asm_sum;

	  	if (tac >= c_asm_sum) {
			nmb5+=0.0001;
			sgn=1;
		} else {
		    nmb5-=0.0001;
			sgn=-1;
		}	 
		   
		if (nmb5<=0) {
			break;
		}
    } while (sgn * presgn > =0);

    for(k = 0; k < age; k++) {
		f[k] = f_asm[k];
	}
}

void output() {
    printf("MED_RPS=%f, MOR=%f\n",MED_RPS,MOR);
    printf("MSY=%f, SSB_MSY=%f, SSB_06=%f, SSB_01=%f\n",MSY,SSB_MSY,SSB_06,SSB_01);
    printf("DELAY=%d, It starts from %d to %d\n",DELAY,START+DELAY,START+DELAY+step);

    switch(type) {
	case 1:
		printf("type 1: RI 1976-1987, RI 1988-2019\n");
        printf("BETA_HIGH=%f, BETA_LOW=%f\n", RI_BETA_HIGH, RI_BETA_LOW);
	    break;
	case 2:
		printf("type 2: BH 1976-1987, BH 1988-2019\n");
        printf("BETA_HIGH=%f, BETA_LOW=%f\n", BH_BETA_HIGH, BH_BETA_LOW);
     	break;
	case 3:
		printf("type 3: RI 1976-1987, HO 1988-1991, RI 1992-2019\n");
        printf("BETA_HIGH=%f, BETA_LOW=%f, BETA_MIDDLE=%f\n", RI_BETA_HIGH_B, HO_BETA_B, RI_BETA_MID_B);
	    break;
	case 4:
		printf("type 4: BH 1976-1987, BH 1988-1991, BH 1992-2019\n");
        printf("BETA_HIGH=%f, BETA_LOW=%f, BETA_MIDDLE=%f\n", BH_BETA_HIGH_B, BH_BETA_LOW_B, BH_BETA_MID_B);
     	break;
    default:
        break;
    }
}

int EstimateFishingMortality() {

	if (SIGMA2 < 0e0) {
		printf("Observation error must be >= 0\n");
		return 1;
	}

	if (SIGMA2 == 0e0) {
		for (k = 0; k < AGE; k++) {
            fishing_mortality_observed[k] = fishing_mortality_true[k];
            abundance[k] = abundance_true[k];
        }
		return 0;
	}

    for (k=0; k < AGE; k++){
        fishing_mortality_observed[k] = fishing_mortality_true[k] * exp(rand_normal(-0.5 * SIGMA2 * SIGMA2, SIGMA2));
        if (fishing_mortality_observed[k] > 3.5e0) {
            fishing_mortality_observed[k]=3.5e0;
        }
		cn[k]=f_Cn(abundance_true[k], MOR, fishing_mortality_true[k]);
        abundance[k]=f_NCn(cn[k], MOR, fishing_mortality_observed[k]);
	}

	return 0;
}

int DetermineCatchWeight() {
	
	if (SIGMA1 < 0e0) {
		return 1;
	}

	for (k = 0; k < AGE; k++){
        if (SIGMA1 > 0e0){   // IE
            fishing_mortality_true[k] = rnd[k] * F_target * exp(rand_normal(-0.5 * SIGMA1 * SIGMA1, SIGMA1));
            if (fishing_mortality_true[k] > 3.5e0) {
                fishing_mortality_true[k] = 3.5e0;
            }
        } else {            // No implemention error
            fishing_mortality_true[k] = rnd[k] * F_target;
        }
  	    catch[k] = f_Cw(abundance_true[k], weight[k], MOR, fishing_mortality_true[k]);
    }
}

void CalculateAbundanceNextyear(double abundance_next_year, double abundance, double natural_mortality, double fishing_mortality) {
    tmp_x1=f_Nl(abundance_true[0], natural_mortality, fishing_mortality[0]);
    tmp_x2=f_Nl(abundance_true[1], natural_mortality, fishing_mortality[1]);
    tmp_x3=f_Nl(abundance_true[2], natural_mortality, fishing_mortality[2]);
    tmp_x4=f_Nl(abundance_true[3], natural_mortality, fishing_mortality[3]);
 	tmp_x5=f_Nh(abundance_true[4], abundance_true[5], natural_mortality, fishing_mortality[4], fishing_mortality[5]);
	abundance_true[1]=tmp_x1;
    abundance_true[2]=tmp_x2;
    abundance_true[3]=tmp_x3;
    abundance_true[4]=tmp_x4;
    abundance_true[5]=tmp_x5;

    ssb_true = 0e0;
	for (k = 1; k < AGE; k++){
	    weight_true[k] = (double)weight_data[ipD+1][k+1];   //  Weight.dat: year 0 1 2 3 4 5
	    ssb_true += f_B(abundance_true[k], mature_rate[k], weight_true[k]);
	}

	fyr=START + ipD + 1;
    rps = CalculateRPS(type, ssb_true, fyr);
	abundance_true[0] = ssb_true * rps;
}

void DetermineTACorTAE() {
 
	if (tac == 0e0){
        for (k = 0; k < AGE; k++) {
	        catch[k] = 0e0;
            fishing_mortality_observed[k] = 0e0;
            fishing_mortality_true[k] = 0e0;
		    if (PHI > 0e0){
		        abundance[k] = abundance_true[k] * exp(rand_normal(-0.5 * PHI * PHI, PHI));
                continue;
		    }
		    abundance[k] = abundance_true[k];
        }
        catch_total = 0e0;
        pattern = 3;
	} else {
	    if ((tac<=tae)||(F_target==3.5e0)){
	        if (SIGMA1>0e0){   // IE
                for (k=0;k<AGE;k++) {
                    rnd[k] = rnd[k] * exp(rand_normal(-0.5 * SIGMA1 * SIGMA1, SIGMA1));
                }
            }
            mult_TAC(tac,rnd,x_t,weigh,f,AGE,MOR); // change fishing_mortality[]
            for (k = 0; k < AGE; k++){
                fishing_mortality_true[k] = fishing_mortality[k];
                catch[k] = f_Cw(abundance_true[k], weight[k], MOR, fishing_mortality_true[k]);
	        }
	        catch_total=tac;
            pattern=2;
	    } else {
	        catch_total=tae;
            pattern=1;
	    }
    
        return_value = EstimateFishingMortality();
        if (return_value != 0) {
            return 1;
        }
    }
}