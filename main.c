#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <fcntl.h>

enum {n = 4};

long double a[n + 1];
long double b[n + 1];
long double sum[n + 1][n + 1];

long double
y_Lagrange(long double *x, long double *y, int num, long double x0) {
    long double L;
    int i, j;
    long double P;

    L = 0;
    for (i = 0; i < num; i++) {
        P = 1;
        for (j = 0; j < num; j++) {
            if (i != j) {
                P *= (x0 - x[j]) / (x[i] - x[j]);              
            }
        }
        L += y[i] * P;
    }
    return L;
}

long double
y_LeastSquares(long double *x, long double* y, int num, long double x0) {
    int i, j, k;

    for (i = 0; i < n + 1; i++){
        a[i]=0;
        b[i]=0;
            for ( j = 0; j < n + 1; j++){
                sum[i][j] = 0;
            }
    }

    for (i = 0; i < n + 1; i++) {
        for (j = 0; j < n + 1; j++) {
            sum[i][j] = 0;
            for (k = 0; k < num; k++) {
                sum[i][j] += pow(x[k], i + j);
            }
        }
    }
    
    for (i = 0; i < n + 1; i++) {
        for (k = 0; k < num; k++) {
            b[i] += pow(x[k], i) * y[k];
        }
    }

    long double temp = 0;
    for (i = 0; i < n + 1; i++) {
        if (sum[i][i] == 0) {
            for (j = 0; j < n + 1; j++) {
                if (j == i) {
                    continue;
                }
                if (sum[j][i] != 0 && sum[i][j] != 0) {
                    for (k = 0; k < n + 1; k++) {
                        temp = sum[j][k];
                        sum[j][k] = sum[i][k];
                        sum[i][k] = temp;
                    }
                    temp = b[j];
                    b[j] = b[i];
                    b[i] = temp;
                    break;
                }
            }
        }
    }

    for (k = 0; k < n + 1; k++) {
        for (i = k + 1; i < n + 1; i++) {
            if (sum[k][k] == 0) {
                printf("\n Solution doesn't exist! \n");
                return 0;
            }
            long double M = sum[i][k] / sum[k][k];
            for (j = k; j < n + 1; j++) {
                sum[i][j] -= M * sum[k][j];
        }
        b[i] -= M * b[k];
    }
}

    for (i = n; i >= 0; i--) {
        long double s = 0;
        for (j = i; j < n + 1; j++){
            s += sum[i][j] * a[j];
        }
        a[i] = (b[i] - s) / sum[i][i];
    }

    long double result = 0;

    for (i = 0; i < n + 1; i++){
        result += a[i] * pow(x0, i);
    }
    return result;
}

int
main(void) {

    FILE* dots = fopen("./dots.txt", "r");
    FILE* input = fopen("./input.txt", "r");
    FILE* out_Lagrange = fopen("./output_Lagrange.txt", "r+");
    FILE* out_LeastSquares = fopen("./output_LeastSquares.txt", "r+");

    int num = 0;
    long double z;
    while (fscanf(dots, "%Lf", &z) > 0) {
        num++;
    }
    num = num / 2;

    long double x[num];
    long double y[num];

    rewind(dots);
    for (int i = 0; i < num; i++) {
        fscanf(dots, "%Lf", &x[i]);
        fscanf(dots, "%Lf", &y[i]);
    }
    fclose(dots);

    long double x0;
    while (fscanf(input, "%Lf", &x0) > 0) {
        fprintf(out_Lagrange, "%Lf %Lf\n", x0, y_Lagrange(x, y, num, x0));
        fprintf(out_LeastSquares, "%Lf %Lf\n", x0, y_LeastSquares(x, y, num, x0));
    }

    printf("Done! You can check output files!\n");

    printf("Least Squares approximate function:\ny = ");
    for (int i = 0; i < n + 1; i++) {
        int sign;
        if (a[i] < 0) {
            sign = -1;
            printf(" - ");
        } else {
            sign = 1;
            if (i != 0) {
                printf(" + ");
            }
        }
        printf("%Lf ", a[i] * sign);
        if (i != 0) {
            printf ("* x^%d", i);
        }
    }
    printf("\n");
    fclose(out_Lagrange);
    fclose(out_LeastSquares);

    return 0;
} 