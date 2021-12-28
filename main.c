#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

enum {n = 4};

double x[11] = {0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2};
double y[11] = {0, 0.006732, 0.058195, 0.030482, 0.387483, 0.958924, 0.48283, 1.802771, 4.052411, 2.403475, 4.352169};
 
double 
y_Lagrange(double x0) {
    double L;
    int i, j;
    double P;

    L = 0;
    for (i = 0; i <= 10; i++) {
        P = 1;
        for (j = 0; j <= 10; j++) {
            if (i != j) {
                P *= (x0 - x[j]) / (x[i] - x[j]); 
                printf("(x - %lf) / (%lf - %lf) *", x[j], x[i], x[j]);               
            }
        }
        printf("%lf + ", y[i]);
        L += y[i] * P;
    }
    return L;
}

double
y_LeastSquares(double x0) {
    double a0, a1, a2, a3, a4;

    return (a0 + a1 * x0 + a2 * pow(x0, 2) + a3 * pow(x0, 3) + a4 * pow(x0, 4));
}

int
main(void) {
    
    double x0;
    printf("Введите значение точки x0: ");
    scanf("%lf", &x0);

    for (int i = 0; i <= 10; i++){
        printf("%lf %lf\n", x[i], y[i]);
    }

    printf("Lagrange: y0 = %lf\n", y_Lagrange(x0));

    return 0;
}