#include <math.h>
#include <stdio.h>

#define max(x,y) ( (x) < (y) ? (y) : (x) )
#define min(x,y) ( (x) < (y) ? (x) : (y) )

#define ATTEMPTS 100
#define MIN_SCALE_FACTOR 0.125
#define MAX_SCALE_FACTOR 4.0
#define _USE_MATH_DEFINES

#define n 12

static double Runge_Kutta(double a, void(*f)(double, double, double*, double*),  double y[][n], double x,
double h);

void f(double a, double x, double *y, double *ans)
{
    ans[0] = y[3];
    ans[1] = y[4];
    ans[2] = y[5];
    
    ans[3] = y[9];
    ans[4] = y[10];
    ans[5] = y[11];
    
    ans[6] = -y[0];
    ans[7] = -y[1];
    ans[8] = -y[2];
    
    ans[9] = -y[6];
    ans[10] = -y[7];
    ans[11] = -y[8];
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
    
    alpha=-y[1][0]/y[1][1] -(1-y[1][3]+(y[1][4]*y[1][0])/y[1][1]) /((y[1][5]*y[1][1]/y[1][2])-y[1][4]);
    beta =y[1][1]* (1-y[1][3]+(y[1][4]*y[1][0])/y[1][1]) /(y[1][5]*y[1][1]-y[1][2]*y[1][4]);
    y[1][0]=y[1][0]+alpha*y[1][1]+beta*y[1][2];
    y[1][1]=y[1][3]+alpha*y[1][4]+beta*y[1][5];
    y[1][2]=y[1][6]+alpha*y[1][7]+beta*y[1][8];
    y[1][3]=y[1][9]+alpha*y[1][10]+beta*y[1][11];
    y[1][4]=alpha;
    y[1][5]=beta;
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
    double h = 0.000001, eps = 1e-9;
    double a = 0.5, alpha,beta;
    double h_next, err;
    double x_start = 0, x_end  = 1;
    double pogr;
    double m= 1,b,c,x=0;
    //printf("%.6f", lambda(a,1));
    y[0][0] = 0;
    y[0][1] = 0;
    y[0][2] = 0;
    
    y[0][3] = 0;
    y[0][4] = 1;
    y[0][5] = 0;
    
    y[0][6] = 0;
    y[0][7] = 0;
    y[0][8] = 1;
    
    y[0][9] = 0;
    y[0][10] = 0;
    y[0][11] = 0;
    
    err = Prince_Dormand(a, alpha, beta, &f, &lambda, y, x_start, h, x_end, &h_next, eps); 
    alpha=y[1][4];
    beta=y[1][5];
    printf("%.15f, %.15f", alpha, beta);

 }

//-0.500000295211506, -300.000391699326542    0.1

//-0.500004761047938, -75.001571207842986 0.2
//-0.502989290781072, -3.039474645049093  1
//-0.551200593688885, -0.920064376872628   2
//138.557303705548094, 131.042445100575378   3.931849978933
//0.815115893599342, 0.794317271575859    5
//-3.389418705820504, -3.389083867852695    10