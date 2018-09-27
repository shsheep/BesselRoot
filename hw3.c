// 2013011060_YANG_SEUNG_HO
/* Driver for routine rtbis */

#include <stdio.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define N 10000 // 100 intervals in [1.0, 10.0]
#define NBMAX 20 // the maximum quantity of roots to find
#define X1 1.0 // range
#define X2 10.0
#define EULER 2.71828
#define PI 3.14159

static float fx(float x)
{
	return bessj0(x);
}

static void funcd(float x,float *fn, float *df)
{
	*fn=bessj0(x);
	*df = -bessj1(x);
}

static float eq_1(float x)
{
	return 10*pow(EULER, -x)*sin(2*PI*x)-2;
}

static void eq_1d(float x, float *fn, float *df)
{
	*fn = 10*pow(EULER, -x)*sin(2*PI*x)-2;
	*df = 10*pow(EULER, -x)*(2*PI*cos(2*PI*x)-sin(2*PI*x));
}

static float eq_2(float x)
{
	return x*x-2*x*pow(EULER, -x)+pow(EULER, -2*x);
}

static void eq_2d(float x, float *fn, float *df)
{
	*fn = x*x-2*x*pow(EULER, -x)+pow(EULER, -2*x);
	*df = 2*(1+pow(EULER, -x))*(x-pow(EULER, -x));
}

static float eq_3(float x)
{
	return cos(x+sqrt(2))+x*(x/2+sqrt(2));
}

static void eq_3d(float x, float *fn, float *df)
{
	*fn = cos(x+sqrt(2))+x*(x/2+sqrt(2));
	*df = -sin(x+sqrt(2))+x+sqrt(2);
}
/*
static float eq_4(float x)
{
	// heart-shaped equation
	return pow((x*x+y*y-1), 3)-x*x*y*y*y;
}
*/
int main(void) {
	int i, nb=NBMAX;
	float xacc, root, *xb1, *xb2;

	// Memory allocation for xb arrays
	// xb1 and xb2 have x value near roots
	xb1=vector(1,NBMAX);
	xb2=vector(1,NBMAX);

	// First use bracketing routine to get the approximate range of each root.
	zbrak(fx,X1,X2,N,xb1,xb2,&nb);

	printf("Finding roots of Bessel function J0\n");
	printf("via range of [%f, %f] with %d intervals", X1, X2, N);
	
	printf("\nBisection method\n%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtbis(fx,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %14.6f\n",i,root,fx(root));
	}

	printf("\nLinear interpolation method\n%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtflsp(fx,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %14.6f\n",i,root,fx(root));
	}

	printf("\nSecant method\n%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtsec(fx,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %14.6f\n",i,root,fx(root));
	}
	
	printf("\nNR with bracketing method\n%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtsafe(funcd,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %14.6f\n",i,root,fx(root));
	}

	printf("\nNewton-Raphson method\n%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtnewt(funcd,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %14.6f\n",i,root,fx(root));
	}

	// Root finding for other equations via different range
	zbrak(eq_1,0.1,1.0,N,xb1,xb2,&nb);
	printf("\nNR with bracketing method 10*pow(EULER, -x)*sin(2*PI*x)-2\n%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtsafe(eq_1d,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %14.6f\n",i,root,fx(root));
	}

	zbrak(eq_2,0.0,1.0,N,xb1,xb2,&nb);
	printf("\nNR with bracketing method x*x-2*x*pow(EULER, -x)+pow(EULER, -2*x)\n%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtsafe(eq_2d,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %14.6f\n",i,root,fx(root));
	}

	zbrak(eq_3,-2.0,-1.0,N,xb1,xb2,&nb);
	printf("\nNR with bracketing method cos(x+sqrt(2))+x*(x/2+sqrt(2))\n%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtsafe(eq_3d,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %14.6f\n",i,root,fx(root));
	}
/*
	zbrak(fx,-2.0,2.0,N,xb1,xb2,&nb);
	printf("\nNR with bracketing method\n%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtsafe(eq_4,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %14.6f\n",i,root,fx(root));
	}
*/
	free_vector(xb2,1,NBMAX);
	free_vector(xb1,1,NBMAX);
}


