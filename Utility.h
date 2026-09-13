#ifndef _UTILITY_H_
#define _UTILITY_H_

int StoreItems(const char *filename, int *array, int rows, int cols);
int StoreItemsDouble(const char *filename, double *array, int rows, int cols);
double CalculateSumOfArray(double *array, int size_array);
void DisplayItems(double *array, int size_array, FILE *fp);

#endif // _UTILITY_H_