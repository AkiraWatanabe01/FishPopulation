// Cohort Analysis
// pCohort_rule_o.c
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include"MT.h"
#include"FishCode.h"
#include"Utility.h"

// =======================
#define MED_RPS (double)(24e0)    // bi/kg
#define SIGMA1 (double)(0e0)
#define SIGMA2 (double)(0e0)
#define PHI (double)(0e0)
#define MOR (double)(0.4)         // mortality
#define F_0 (double)(0.18)
#define F_1 (double)(0.18)
#define F_2 (double)(0.24)
#define F_3 (double)(0.5)
#define F_4 (double)(0.5)
#define F_5 (double)(0.5)
#define RHO_0 (double)(0.149)
#define RHO_1 (double)(0.407)
#define RHO_2 (double)(0.435)
#define RHO_3 (double)(0.551)
#define RHO_4 (double)(1e0)
#define RHO_5 (double)(1e0)
#define SSB_MSY (double)(1187e0)
#define SSB_06 (double)(487e0)
#define SSB_01 (double)(69e0)
#define MSY (double)(389e0)       // kilo ton
#define YEAR 43
#define START 1976
#define TRIAL 1000
#define AGE 6

int main(int argc, char *argv[]){
    time_t start_time, end_time;
    double abundance_first_year[AGE];
    double abundance_second_year[AGE];
    double catch_first_year[AGE];
    double catch_second_year[AGE];
    double fishing_mortality_first_year[AGE];
    double fishing_mortality_second_year[AGE];
    double tmp_x0, tmp_x1, tmp_x2, tmp_x3, tmp_x4, tmp_x5;
    double abundance_total;
    double biomass_total;
    double catch_total;
    double biomass_true_total;
    double ssb, rps;
    double coef;
    double fishing_mortality_observed[AGE];
    double fishing_mortality_previous_year[5][AGE];   // fishing_mortality_previous_year[year 0:4][age 0:5]
    double abc1, abc2;
    double fishing_mortality_true[AGE];
    double abundance_true[AGE];
    double biomass_true[AGE];
    double ssb_true;                                  // true number (million)
    double tae, tac, max_rnd;
    double abundance[AGE];
    double fishing_mortality[AGE];
    double catch[AGE];
    double cn[AGE];
    double biomass[AGE];
    double weight[AGE], t1_x[AGE], t2_x[AGE], weight_true[AGE];
    double rnd[AGE];
    double tar_sum, c_asm_sum, presum;
    double F_target;
    double alp[YEAR+1];
    double abs_c, pre_csum, sum_abs_c, sum_deno_c;                     // === AAV ===
    double catch_data[YEAR + 1][7];
    double fishing_mortality_data[YEAR + 1][8];
    double abundance_data[YEAR + 1][7];
    const double mature_rate[AGE] = {0e0, 2e-1, 1e0, 1e0, 1e0, 1e0};
    int DELAY;
    int weight_data[YEAR+1][7];                                        // weight(g)
    int i, j, k, l;
    int ipD;                                                           // i + DELAY
    int step, type;
    int return_value;
    int size_of_column;
    int size_of_row;
    short fyr, pattern;
    char filename[20];
    FILE *fp;

    start_time = time(NULL);
    init_genrand(0);
    DELAY = atoi(argv[1]) - 1976;
    if (DELAY < 0) {
        printf("The first argument is invalid!\n");
        return 1;
    }
    F_target = atof(argv[2]);
    if (F_target <= 0e0) {
        printf("The second argument is invalid!\n");
        return 1;
    }
    type = atoi(argv[3]);

// ============ Abundance =============
    return_value = StoreItemsDouble("../../Data/Abundance.dat", &abundance_data[0][0], (YEAR + 1), 7);
    if (return_value != 0) {
        printf("Error termination\n");
        return 1;
    }

    for (i = 0; i < AGE; i++) {
        abundance_first_year[i] = abundance_data[DELAY][i + 1];
        abundance_second_year[i] = abundance_data[DELAY+1][i + 1];
    }
// ============ end abundance ==============
// ============ Catch ============
    return_value = StoreItemsDouble("../../Data/Catch.dat", &catch_data[0][0], (YEAR + 1), 7);
    if (return_value != 0) {
        printf("Error termination\n");
        return 1;
    }

    for (i = 0; i < AGE; i++) {
        catch_first_year[i] = catch_data[DELAY][i + 1];
        catch_second_year[i] = catch_data[DELAY+1][i + 1];
    }
//    ============ end Catch ==============
//    ============ Fishing mortality =============
    return_value = StoreItemsDouble("../../Data/FishingMortality.dat", &fishing_mortality_data[0][0], (YEAR + 1), 8);
    if (return_value != 0) {
        printf("Error termination\n");
        return 1;
    }

    for (i = 0; i < AGE; i++) {
        fishing_mortality_first_year[i] = fishing_mortality_data[DELAY][i + 1];
        fishing_mortality_second_year[i] = fishing_mortality_data[DELAY+1][i + 1];
    }
//     ============ end fishing mortality ===============
//     ============ Alpha ================
    size_of_row = sizeof(alp) / sizeof(alp[0]);
    switch(type) {
    case 2:
        return_value = StoreItemsDouble("../../Data/Alpha.dat", &alp[0], size_of_row, 1);
        if (return_value != 0) {
            printf("Error termination\n");
            return 1;
        }
        break;
    case 4:
        return_value = StoreItemsDouble("../../Data/Alpha3.dat", &alp[0], size_of_row, 1);
        if (return_value != 0) {
            printf("Error termination\n");
            return 1;
        }
        break;
    default:
        break;
    }
// ============ end Alpha =============
// ============ Weight =============
    return_value = StoreItems("../../Data/Weight.dat", &weight_data[0][0], (YEAR + 1), 7);
    if (return_value != 0) {
        printf("Error termination\n");
        return 1;
    }
// ============ end Weight =============

    fp = fopen("Cohort_rule.dat", "w");
// ============= begin TRIAL ==================
    abc2 = 0e0;
    for (i = 0; i < AGE; i++) {
        fishing_mortality[i] = fishing_mortality_first_year[i];
        abundance[i] = abundance_first_year[i];
        catch[i] = catch_first_year[i];
    }

// ============
    for(l = 0; l < 5; l++){
        fishing_mortality_previous_year[l][0] = RHO_0;
        fishing_mortality_previous_year[l][1] = RHO_1;
        fishing_mortality_previous_year[l][2] = RHO_2;
        fishing_mortality_previous_year[l][3] = RHO_3;
        fishing_mortality_previous_year[l][4] = RHO_4;
        fishing_mortality_previous_year[l][5] = RHO_5;
    }

    sum_abs_c = 0e0;
    sum_deno_c = 0e0;
// ============= begin 1st and 2nd year ===================
    for (i = 0; i <= 1; i++){
        ipD = i + DELAY;

        ssb = 0e0;

        for (k = 0; k < AGE; k++){
    	    weight[k] = (double)weight_data[ipD][k + 1];
    	    biomass[k] = f_W(abundance[k], weight[k]);
            ssb += f_B(abundance[k], mature_rate[k], weight[k]);

            fishing_mortality_previous_year[4][k] = fishing_mortality_previous_year[3][k];
    	    fishing_mortality_previous_year[3][k] = fishing_mortality_previous_year[2][k];
    	    fishing_mortality_previous_year[2][k] = fishing_mortality_previous_year[1][k];
    	    fishing_mortality_previous_year[1][k] = fishing_mortality_previous_year[0][k];
            fishing_mortality_previous_year[0][k] = fishing_mortality[k];
        }

	    biomass_total = CalculateSumOfArray(&biomass[0], (sizeof(biomass) / sizeof(biomass[0])));
        abundance_total = CalculateSumOfArray(&abundance[0], (sizeof(abundance) / sizeof(abundance[0])));
        catch_total = CalculateSumOfArray(&catch[0], (sizeof(catch) / sizeof(catch[0])));

        printf(" %4d",START+ipD);
        DisplayItems(&fishing_mortality[0], (sizeof(fishing_mortality) / sizeof(fishing_mortality[0])), stdout);
        DisplayItems(&biomass[0], (sizeof(biomass) / sizeof(biomass[0])), stdout);
        printf(" %.2f",biomass_total);
        DisplayItems(&catch[0], (sizeof(catch) / sizeof(catch[0])), stdout);
        printf(" %.2f",catch_total);
        printf("| %f %f %f", coef, ssb, abundance[0] / ssb);
        printf("| %f %f %d", 0e0, 0e0, 0);
        printf("| %d", 0);
        printf("\n");

        fprintf(fp," %4d", START + i + DELAY);
        DisplayItems(&fishing_mortality[0], (sizeof(fishing_mortality) / sizeof(fishing_mortality[0])), fp);
        DisplayItems(&biomass[0], (sizeof(biomass) / sizeof(biomass[0])), fp);
        fprintf(fp," %.2f", biomass_total);
        DisplayItems(&catch[0], (sizeof(catch) / sizeof(catch[0])), fp);
        fprintf(fp," %.2f",catch_total);
        fprintf(fp," %f %f %f",coef,ssb,abundance[0]/ssb);
        fprintf(fp," %f %f %d",0e0,0e0,0);
        fprintf(fp," %d",0);
        fprintf(fp,"\n");

// ======== begin t+1 =========
        CalculateAbundanceNextyear(&t1_x[0], &abundance[0], MOR, &fishing_mortality[0]);
        //ssb=0e0;
        //for (k = 1; k < AGE; k++) {
        //    ssb+=f_B(t1_x[k], mature_rate[k], weight[k]);
        //}
    	//t1_x[0] = ssb * MED_RPS;

// ====== end t+1 ======
// ====== start t+2 ======
// ============ begin f_current =============
// !!!!! fishing_mortality[year 0:4][age 0:5]  !!!!!
    	for (k = 0; k < AGE; k++) {
            fishing_mortality[k] = 0.2 * (fishing_mortality_previous_year[0][k] + fishing_mortality_previous_year[1][k] + fishing_mortality_previous_year[2][k] + fishing_mortality_previous_year[3][k] + fishing_mortality_previous_year[4][k]);
        }
// ============ end f_current =============
        CalculateAbundanceNextyear(&t2_x[0], &t1_x[0], MOR, &fishing_mortality[0]);
        //ssb=0e0;
	    //for (k = 1; k < AGE; k++) {
        //    ssb += f_B(t2_x[k], mature_rate[k], weight[k]);
        //}
	    //t2_x[0]=ssb*MED_RPS;

        coef = CalculateCoefficientForABC(ssb, SSB_01, SSB_MSY, SSB_06);
        fishing_mortality[0] = coef * F_0;
        fishing_mortality[1] = coef * F_1;
        fishing_mortality[2] = coef * F_2;
        fishing_mortality[3] = coef * F_3;
        fishing_mortality[4] = coef * F_4;
        fishing_mortality[5] = coef * F_5;

	    tar_sum=0e0;
        for (k = 0; k < AGE; k++) {
            tar_sum+=f_Cw(t2_x[k], weight[k], MOR, fishing_mortality[k]);
        }
// ========== end t+2 ============

        abc1 = abc2;
        abc2 = tar_sum;

        for (j = 0; j < AGE; j++) {
            fishing_mortality[j] = fishing_mortality_second_year[j];
        	catch[j] = catch_second_year[j];
        	abundance[j] = abundance_second_year[j];
        }

    	printf("ABC_t+2:");
    	printf(" %.6f",abc2);
    	printf("\n");
    }
// ============= end 1st and 2nd year ===================

// =============== true x at 2 ============
    ssb_true = 0e0;
    for (k = 1; k < AGE; k++){
        abundance_true[k] = t1_x[k];
        weight_true[k] = (double)weight_data[DELAY+2][k + 1];   //  Weight.dat: year 0 1 2 3 4 5
        ssb_true += f_B(abundance_true[k], mature_rate[k], weight_true[k]);
    }

    fyr = START + DELAY + 2;
    rps = CalculateRPS(type, ssb_true, fyr);
    abundance_true[0] = ssb_true * rps;
    pre_csum = catch_total;
// ===============              =============
    printf("\n");
    step = YEAR - DELAY;
    for (i = 2; i <= step; i++){
        ipD = i + DELAY;
        printf("START Year %4d, ABC:%f\n", START + ipD, abc1);
        printf("rps=%f, ssb_true=%f, abundance[0]=%f, abundance_true[0]=%f\n", rps, ssb_true, abundance[0], abundance_true[0]);
        // =========== begin t ===========
        for (k = 0; k < AGE; k++) {
            weight[k] = (double)weight_data[ipD][k+1];   //  Weight.dat: year 0 1 2 3 4 5
            rnd[0] = RHO_0;
            rnd[1] = RHO_1;
            rnd[2] = RHO_2;
	        rnd[3] = RHO_3;
            rnd[4] = RHO_3;
            rnd[5] = RHO_3;
        }

        return_value = DetermineCatchWeight();
        if (return_value != 0) {
            return 1;
        }
        tae = CalculateSumOfArray(&catch[0], (sizeof(catch) / sizeof(catch[0])));

        tac=0e0;
        for (k = 0; k < AGE; k++) {
            tac += f_Cw(abundance_true[k], weight[k], MOR, 3.5e0);
        }
        if (abc1 <= tac) {
            tac = abc1;
        }
        printf("TAC: %.2f\n", tac);
        printf("TAE: %.2f\n", tae);

        // ==================== min{tac, tae} ===============================
        DetermineTACorTAE();
    }

    biomass_true_total = 0e0;
    abundance_total = 0e0;
    ssb_true = 0e0;
 	for (k = 0; k < AGE; k++) {
	    biomass_true[k] = f_W(abundance_true[k], weight[k]);
  	    biomass_true_total += biomass_true[k];
        abundance_total += abundance_true[k];
        ssb_true += f_B(abundance_true[k], mature_rate[k], weight_true[k]);
        // !!!!! fishing_mortality[year 0:4][age 0:5]  !!!!!
        fishing_mortality_previous_year[4][k] = fishing_mortality_previous_year[3][k];
        fishing_mortality_previous_year[3][k] = fishing_mortality_previous_year[2][k];
        fishing_mortality_previous_year[2][k] = fishing_mortality_previous_year[1][k];
        fishing_mortality_previous_year[1][k] = fishing_mortality_previous_year[0][k];
        fishing_mortality_previous_year[0][k] = fishing_mortality_observed[k];
	}

// =========== begin AAV ===========
    abs_c = catch_total - pre_csum;
	if (abs_c<0) abs_c=-1e0*abs_c;
	sum_abs_c+=abs_c;                // numerator
	sum_deno_c+=pre_csum;
    printf("AAV: %.2f %.2f %.2f %.2f %.4f\n", catch_total, pre_csum, sum_abs_c, sum_deno_c, sum_abs_c / sum_deno_c);
    pre_csum = catch_total;
// =========== end AAV ===========
///*
    printf(" %4d",START+i+DELAY);
    DisplayItems(&fishing_mortality_true[0], (sizeof(fishing_mortality_true) / sizeof(fishing_mortality_true[0])), stdout);
    DisplayItems(&biomass_true[0], (sizeof(biomass_true) / sizeof(biomass_true[0])), stdout);
    printf(" %.2f",biomass_true_total);
    DisplayItems(&catch[0], (sizeof(catch) / sizeof(catch[0])), stdout);
    printf(" %.2f",catch_total);
    printf("| %f %f %f",coef,ssb_true,rps);
    printf(" %f %f %d",tac,tae,pattern);
	printf(" %.4f",sum_abs_c/sum_deno_c);
    printf("\n");
//*/

  	fprintf(fp," %4d",START+i+DELAY);
    DisplayItems(&fishing_mortality_true[0], (sizeof(fishing_mortality_true) / sizeof(fishing_mortality_true[0])), fp);
    DisplayItems(&biomass_true[0], (sizeof(biomass_true) / sizeof(biomass_true[0])), fp);
  	fprintf(fp," %.2f",biomass_true_total);
    DisplayItems(&catch[0], (sizeof(catch) / sizeof(catch[0])), fp);
    fprintf(fp," %.2f",catch_total);
  	fprintf(fp," %f %f %f",coef,ssb_true,rps);
  	fprintf(fp," %f %f %d",tac,tae,pattern);
	fprintf(fp," %.4f",sum_abs_c/sum_deno_c);
  	fprintf(fp,"\n");
// =========== end t ===========
// ======== begin t+1 =========
    CalculateAbundanceNextyear(&t1_x[0], &abundance[0], MOR, &fishing_mortality_observed[0]);

// ==== variable weight end  ====
//    ssb=0e0;
//	for (k=1;k<AGE;k++)  ssb += f_B(t1_x[k], mature_rate[k], weight[k]);
//	t1_x[0]=ssb * MED_RPS;

// =============== true x at t+1 ============
    CalculateAbundanceNextyear(abundance_next, &abundance_true[0], MOR, &fishing_mortality_true[0]);
// ======== end t+1 =========
// ======== begin t+2 =========
    // ============ begin f_current =============
    // !!!!! fishing_mortality[year 0:4][age 0:5]  !!!!!
	for (k = 0; k < AGE; k++) {
        fishing_mortality[k] = 0.2 * (fishing_mortality_previous_year[0][k] + fishing_mortality_previous_year[1][k] + \
                               fishing_mortality_previous_year[2][k] + fishing_mortality_previous_year[3][k] + \
                               fishing_mortality_previous_year[4][k]);
    }
    // ============ end f_current =============
    CalculateAbundanceNextyear(&t2_x[0], &t1_x[0], MOR, &fishing_mortality_true[0]);

	//ssb = 0e0;
	//for (k = 1; k < AGE; k++) {
    //    ssb += f_B(t2_x[k], mature_rate[k], weight[k]);
    //}
	//t2_x[0] = ssb * MED_RPS;

    coef = CalculateCoefficientForABC(ssb, SSB_01, SSB_MSY, SSB_06);
	fishing_mortality[0] = coef * F_0;
    fishing_mortality[1] = coef * F_1;
    fishing_mortality[2] = coef * F_2;
    fishing_mortality[3] = coef * F_3;
    fishing_mortality[4] = coef * F_4;
    fishing_mortality[5] = coef * F_5;
// ======== end t+2 =========
// ============ begin ABC t+2 ==============
    tar_sum = 0e0;
    for (k = 0; k < AGE; k++) {
        tar_sum += f_Cw(t2_x[k], weight[k], MOR, fishing_mortality[k]);
    }
    abc1 = abc2;
    abc2 = tar_sum;
    printf("ABC_t+2: %.2f\n", abc2);
    printf("\n");
// ============ end ABC t+2 ==============
  }
// ============ end TRIAL ====================
    fclose(fp);
    end_time = time(NULL);

    printf("time:%ld\n", end_time - start_time);

    return 0;
}
