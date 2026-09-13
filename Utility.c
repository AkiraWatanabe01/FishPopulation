#include <stdio.h>
#include <stdlib.h>

int StoreItems(const char *filename, int *array, int rows, int cols) {
    FILE *fp;
    int i, j;
    int return_fscanf;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("The file does not open: %s\n", filename);
        return 1;
    }

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            return_fscanf = fscanf(fp, "%d", &array[i * cols + j]);
            if (return_fscanf != 1) {
                printf("Error occures when reading data.\n");
                fclose(fp);
                return 1;
            }
        }
    }

    fclose(fp);
    return 0;
}

int StoreItemsDouble(const char *filename, double *array, int rows, int cols) {
    FILE *fp;
    int i, j;
    int return_fscanf;

    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("The file does not open: %s\n", filename);
        return 1;
    }

    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            return_fscanf = fscanf(fp, "%lf", &array[i * cols + j]);
            if (return_fscanf != 1) {
                printf("Error occures when reading data.\n");
                fclose(fp);
                return 1;
            }
        }
    }

    fclose(fp);
    return 0;
}

double CalculateSumOfArray(double *array, int size_array) {
    int i;
    double total = 0e0;

    for (i = 0; i < size_array; i++) {
        total += array[i];
    }

    return total;
}

void DisplayItems(double *array, int size_array, FILE *fp) {
    int k;

    for (k = 0; k < size_array; k++) {
        fprintf(fp," %.2f", array[k]);
    }
}

double rand_normal( double mu, double sigma ){  //gauss random number
    return (mu + sigma * sqrt(-2.0 * log(genrand_real3())) * sin(2.0 * M_PI * genrand_real3()));
}