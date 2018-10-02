
#include <math.h>
#define MAXIT 30

// functions that set Muller constants used in the Muller Algorithm
float set_muller_a (float (*func)(float), float p0, float p1, float p2) {
	float fp0 = (*func)(p0), fp1 = (*func)(p1), fp2 = (*func)(p2);
	float a = ((p1-p2)*(fp0-fp2)-(p0-p2)*(fp1-fp2)) / (p0-p2)*(p1-p2)*(p0-p1);
	return a;
}

float set_muller_b (float (*func)(float), float p0, float p1, float p2) {
	float fp0 = (*func)(p0), fp1 = (*func)(p1), fp2 = (*func)(p2);
	float b = ((p0-p2)*(p0-p2)*(fp1-fp2)-(p1-p2)*(p1-p2)*(fp0-fp2)) / (p0-p2)*(p1-p2)*(p0-p1);
	return b;
}

float rtmuller(float (*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float fl,fh,xl,xh,f,rtf, dx, xm, fm, a, b, c;

	xl = x1;
	xm = x1+2.0;//(x1+x2)/2.0; // set initial third point
	xh = x2;
	fl=(*func)(x1);
	fm = (*func)(xm);
	fh=(*func)(x2);

	for (j=1;j<=MAXIT;j++) {
		// Set Muller constants
		c = fh;
		b = set_muller_b((*func), xl, xm, xh);
		a = set_muller_a((*func), xl, xm, xh);
		// Evaluate the Pn using Pn-1, Pn-2, Pn-3
		dx = (b > 0) ? -2*c/(b+sqrt(b*b-4*a*c)) : ( (b == 0 ) ? -2*c/(sqrt(-4*a*c)) : -2*c/(b-sqrt(b*b-4*a*c)) );
		rtf = xh + dx;

		// For Debug : printf("Before -> a : %e b : %e c : %e xl : %f xm : %f xh : %f dx : %f\n", a, b, c, xl, xm, xh, dx);
		// only update f value to compare it with 0
		f = (*func)(rtf);
		xl = xm;
		xm = xh;
		xh = rtf;

		// For Debug : printf("After -> xl : %f xm : %f xh : %f\n", xl, xm, xh);
		if (fabs(dx) < xacc || f == 0.0) return rtf;
	}
	nrerror("Maximum number of iterations exceeded in rtmuller");
	return 0.0;
}
#undef MAXIT
