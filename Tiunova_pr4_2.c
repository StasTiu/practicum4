#include <math.h>
#include <stdio.h>

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

#define ATTEMPTS 10
#define MIN_SCALE_FACTOR 0.125
#define MAX_SCALE_FACTOR 4.0
#define _USE_MATH_DEFINES

#define n 5

static double Runge_Kutta(double a, void(*f)(double, double, double*, double*),  double y[][n], double x,
double h);

void f(double a, double x, double *y, double *ans)
{
    ans[0] = y[1];
    ans[1] = y[3];
    ans[2] = -y[0];
    ans[3] = -y[2];
    ans[4] = y[3]*y[3]-y[0]*y[0];
    
}
double lambda(double a, double x)
{

    return 2;
}
double Prince_Dormand(double a, double alpha, double beta, void(*f)(double, double, double*,double*),double(*lambda)(double, double x),double y[][n]
,double x, double h, double xmax, double *h_next, double tolerance) {
    double scale,integ[6],integr;
    double temp_y[2][n];
    double err = 0, global_err = 0;
    double l0 = 0, l1 = 0;
    double yy = 0;
    int i, j;
    alpha=2;
    int last_interval = 0;
    if (xmax < x || h <= 0.0) return -2;
    *h_next = h;
    for(i = 0; i < n; i++)
        y[1][i] = y[0][i];
    if (xmax == x) return 0;
    h = min(h, xmax - x);
    tolerance /= (xmax - x);
    for (i = 0; i < n; i++)
        temp_y[0][i] = y[0][i];
    while (x < xmax) {
        scale = 1.0;
        for (i = 0; i < ATTEMPTS; i++) {
            yy = 0;
            err = fabs(Runge_Kutta(a, f, temp_y, x, h));
            if (err == 0.0) { scale = MAX_SCALE_FACTOR; break; }
            for (j = 0; j < n; j++) yy += (temp_y[0][j] == 0.0) ? tolerance:
            fabs(temp_y[0][j]);
            scale = 0.8 * sqrt(sqrt(tolerance * yy / err));
            scale = min(max(scale, MIN_SCALE_FACTOR), MAX_SCALE_FACTOR);
            if (err < (tolerance * yy)) break;
            h *= scale;
            if (x + h > xmax) h = xmax - x;
            else if (x + h + 0.5 * h > xmax) h = 0.5 * h;
        }
        
        if (i >= ATTEMPTS) { *h_next = h * scale; return -1; };
        for (j = 0; j < n; j++) temp_y[0][j] = temp_y[1][j];
        x += h;
        l1 = lambda(a, x);
        global_err = err + global_err * h * max(l0, l1);
        l0 = l1;
        h *= scale;
        *h_next = h;
        if (last_interval) break;
        if (x + h > xmax) { last_interval = 1; h = xmax - x; }
        else if (x + h + 0.5 * h > xmax) h = 0.5 * h;
    }
    for (j = 0; j < n; j++) y[1][j] = temp_y[1][j];
    

    //printf("%.6f, %.6f \n ", alpha, beta);

    
    
    
    return global_err;
}


static double Runge_Kutta(double a, void(*f)(double, double, double*, double*)
, double y[][n], double x0, double h) {
    static const double r_45 = 1.0 / 45.0;
    static const double r_8_9 = 8.0 / 9.0;
    static const double r_6561 = 1.0 / 6561.0;
    static const double r_167904 = 1.0 / 167904.0;
    static const double r_142464 = 1.0 / 142464.0;
    static const double r_21369600 = 1.0 / 21369600.0;
    double y_tmp[n];
    double err = 0;
    double k1[n], k2[n], k3[n], k4[n], k5[n], k6[n], k7[n];
    double h5 = 0.2 * h;
    for (int i = 0; i < n; i++)
        y_tmp[i] = y[0][i];
    (*f)(a, x0, y_tmp, k1);
    for (int i = 0; i < n; i++)
        y_tmp[i] = y[0][i] + h5 * k1[i];
	(*f)(a, x0 + h5, y_tmp, k2);
    for (int i = 0; i < n; i++)
        y_tmp[i] = y[0][i] + h * (0.075 * k1[i] + 0.225 * k2[i]);
    (*f)(a, x0 + 0.3*h, y_tmp, k3);
    for (int i = 0; i < n; i++)
        y_tmp[i] = y[0][i] + h * r_45 * (44.0 * k1[i] - 168.0 * k2[i] + 160 * k3[i]);
    (*f)(a, x0 + 0.8*h, y_tmp, k4);
    for (int i = 0; i < n; i++)
        y_tmp[i] = y[0][i] + r_6561 * h * (19372.0 * k1[i]
        - 76080.0 * k2[i] + 64448.0 * k3[i] - 1908.0 * k4[i]);
    (*f)(a, x0 + r_8_9 * h, y_tmp, k5);
    for (int i = 0; i < n; i++)
        y_tmp[i] = y[0][i] + r_167904 * h * (477901.0 * k1[i] - 1806240.0 * k2[i]
        + 1495424.0 * k3[i] + 46746.0 * k4[i] - 45927.0 * k5[i]);
    (*f)(a, x0 + h, y_tmp, k6);
    for (int i = 0; i < n; i++)
        y_tmp[i] = y[0][i] + r_142464 * h * (12985.0 * k1[i] + 64000.0 * k3[i]
        + 92750.0 * k4[i] - 45927.0 * k5[i] + 18656.0 * k6[i]);
    (*f)(a, x0 + h, y_tmp, k7);
    for (int i = 0; i < n; i++)
    {
        y[1][i] = y[0][i] + r_21369600 * h * (1921409.0 * k1[i] + 9690880.0 * k3[i]
        + 13122270.0 * k4[i] - 5802111.0 * k5[i] + 1902912.0 * k6[i] + 534240.0 * k7[i]);
        err += fabs(r_21369600 * (26341.0 * k1[i] - 90880.0 * k3[i] + 790230.0 * k4[i]
        - 1086939.0 * k5[i] + 895488.0 * k6[i] - 534240.0 * k7[i]));
    }
    return err;
}

int main()
{
    double y[2][n];
    double y1[2][n];
	double y2[2][n];
	double err1, err2;
    double h = 0.000000085, eps;
    double a = 0.5, alpha,beta;
    double h_next, err;
    double x_start = 0, x_end  = 10;
    double pogr;
    double m= 1,b,c,x=0;
    //printf("%.6f", lambda(a,1));
    
    
    eps =1e-5;

    x_end=1;
    
    y[0][0] = 0;
    y[0][1] = -0.502989290781072;
    y[0][2] = -3.039474645049093;
    y[0][3] = 0;
    y[0][4] = 0;
    
    y1[0][0] = 0;
    y1[0][1] = -0.502989290781072;
    y1[0][2] =-3.039474645049093;
    y1[0][3] = 0;
    y1[0][4] = 0;
    
    y2[0][0] = 0;
    y2[0][1] =-0.502989290781072;
    y2[0][2] = -3.039474645049093;
    y2[0][3] = 0;
    y2[0][4] = 0;

    err = Prince_Dormand(a, alpha, beta, &f, &lambda, y, x_start, h, x_end, &h_next, eps); 
    err1 = Prince_Dormand(a, alpha, beta,  &f, &lambda, y1, x_start, h, x_end, &h_next, eps*1e-2);
    err2 = Prince_Dormand(a, alpha, beta, &f, &lambda, y2, x_start, h, x_end, &h_next, eps*1e-4);
    printf("$T=1$ & %.9f & %.9f & %.9f ", y[0][1], y[0][2], y[1][4]);
    printf("& %.9f ", y1[1][4]);
     printf("& %.9f ", y2[1][4]);
    printf("\\\\ \n");
     // printf("$T=1$ & %.9f & %.9f & %.9f \\\\ \n ", y[1][1]-y1[1][1], y1[1][1]-y2[1][1], (y[1][1]-y1[1][1])/(y1[1][1]-y2[1][1]));


    x_end=2;
    y[0][0] = 0;
    y[0][1] = -0.551200593688885;
    y[0][2] = -0.920064376872628 ;
    y[0][3] = 0;
    y[0][4] = 0;

    y1[0][0] = 0;
    y1[0][1] =-0.551200593688885;
    y1[0][2] = -0.920064376872628;
    y1[0][3] = 0;
    y1[0][4] = 0;

    y2[0][0] = 0;
    y2[0][1] = -0.551200593688885;
    y2[0][2] =  -0.920064376872628;
    y2[0][3] = 0;
    y2[0][4] = 0;

    err = Prince_Dormand(a, alpha, beta, &f, &lambda, y, x_start, h, x_end, &h_next, eps);
    err1 = Prince_Dormand(a, alpha, beta,  &f, &lambda, y1, x_start, h, x_end, &h_next, eps*1e-2);
    err2 = Prince_Dormand(a, alpha, beta, &f, &lambda, y2, x_start, h, x_end, &h_next, eps*1e-4);

    printf("$T=2$ & %.9f & %.9f & %.9f ", y[0][1], y[0][2], y[1][4]);
    printf("& %.9f ", y1[1][4]);
    printf("& %.9f ", y2[1][4]);
    printf("\\\\ \n");
    //printf("$T=2$ & %.9f & %.9f & %.9f \\\\ \n", y[1][1]-y1[1][1], y1[1][1]-y2[1][1], (y[1][1]-y1[1][1])/(y1[1][1]-y2[1][1]));



    x_end= 3.931849978933;
    y[0][0] = 0;
    y[0][1] = 138.557303705548094;
    y[0][2] = 131.042445100575378;
    y[0][3] = 0;
    y[0][4] = 0;

    y1[0][0] = 0;
    y1[0][1] = 138.557303705548094;
    y1[0][2] = 131.042445100575378;
    y1[0][3] = 0;
    y1[0][4] = 0;

    y2[0][0] = 0;
    y2[0][1] =138.557303705548094;
    y2[0][2] = 131.042445100575378;
    y2[0][3] = 0;
    y2[0][4] = 0;
    err = Prince_Dormand(a, alpha, beta, &f, &lambda, y, x_start, h, x_end, &h_next, eps);
    err1 = Prince_Dormand(a, alpha, beta,  &f, &lambda, y1, x_start, h, x_end, &h_next, eps*1e-2);
    err2 = Prince_Dormand(a, alpha, beta, &f, &lambda, y2, x_start, h, x_end, &h_next, eps*1e-4);

    printf("$T= 3.931849978933$ & %.9f & %.9f & %.9f ", y[0][1], y[0][2], y[1][4]);
    printf("& %.9f ", y1[1][4]);
    printf("& %.9f ", y2[1][4]);
    printf("\\\\ \n");
    //printf("$T= 3.931849978933$ & %.9f & %.9f & %.9f \\\\ \n", y[1][1]-y1[1][1], y1[1][1]-y2[1][1], (y[1][1]-y1[1][1])/(y1[1][1]-y2[1][1]));

    x_end=5;
    y[0][0] = 0;
    y[0][1] = 0.815115893599342;
    y[0][2] =  0.794317271575859;
    y[0][3] = 0;
    y[0][4] = 0;

    y1[0][0] = 0;
    y1[0][1] = 0.815115893599342;
    y1[0][2] = 0.794317271575859 ;
    y1[0][3] = 0;
    y1[0][4] = 0;

    y2[0][0] = 0;
    y2[0][1] = 0.815115893599342;
    y2[0][2] = 0.794317271575859 ;
    y2[0][3] = 0;
    y2[0][4] = 0;

    err = Prince_Dormand(a, alpha, beta, &f, &lambda, y, x_start, h, x_end, &h_next, eps);
    err1 = Prince_Dormand(a, alpha, beta,  &f, &lambda, y1, x_start, h, x_end, &h_next, eps*1e-2);
    err2 = Prince_Dormand(a, alpha, beta, &f, &lambda, y2, x_start, h, x_end, &h_next, eps*1e-4);

    printf("$T=5$ & %.9f & %.9f & %.9f ", y[0][1], y[0][2], y[1][4]);
    printf("& %.9f ", y1[1][4]);
    printf("& %.9f ", y2[1][4]);
    printf("\\\\ \n");
    //printf("$T=5$ & %.9f & %.9f & %.9f \\\\ \n", y[1][1]-y1[1][1], y1[1][1]-y2[1][1], (y[1][1]-y1[1][1])/(y1[1][1]-y2[1][1]));


    x_end=10;
    y[0][0] = 0;
    y[0][1] =-3.389418705820504;
    y[0][2] = -3.389083867852695;
    y[0][3] = 0;
    y[0][4] = 0;
    
    y1[0][0] = 0;
    y1[0][1] = -3.389418705820504;
    y1[0][2] = -3.389083867852695;
    y1[0][3] = 0;
    y1[0][4] = 0;
    
    y2[0][0] = 0;
    y2[0][1] = -3.389418705820504;
    y2[0][2] = -3.389083867852695;
    y2[0][3] = 0;
    y2[0][4] = 0;
    err = Prince_Dormand(a, alpha, beta, &f, &lambda, y, x_start, h, x_end, &h_next, eps); 
    err1 = Prince_Dormand(a, alpha, beta,  &f, &lambda, y1, x_start, h, x_end, &h_next, eps*1e-2);
    err2 = Prince_Dormand(a, alpha, beta, &f, &lambda, y2, x_start, h, x_end, &h_next, eps*1e-4);
    
    printf("$T=10$ & %.9f & %.9f & %.9f ", y[0][1], y[0][2], y[1][4]);
    printf("& %.9f ", y1[1][4]);
    printf("& %.9f ", y2[1][4]);
    printf("\\\\ \n");
    //printf("$T=10$ & %.9f & %.9f & %.9f \\\\ \n ", y[1][1]-y1[1][1], y1[1][1]-y2[1][1], (y[1][1]-y1[1][1])/(y1[1][1]-y2[1][1]));


 }

