
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
	float fl,fh,xl,xh,swap,del,f,rtf, xm, fm, a, b, c; //dx

	xm = (x1+x2)/2.0;
	fl=(*func)(x1);
	fm = (*func)(xm);
	fh=(*func)(x2);


	// if ( fh )
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xl=x2;
		`xh=x1;
		swap=fl;
		fl=fh;
		fh=swap;
	}
	//dx=xh-xl;
	for (j=1;j<=MAXIT;j++) {
		c = fh;
		b = set_muller_b((*func), xl, xm, xh);
		a = set_muller_a((*func), xl, xm, xh);
		rtf = xh - 2*c/(b+sgn(b)*sqrt(b*b-4*a*c));
		f=(*func)(rtf);
		if (f < 0.0) {
			del=xl-rtf;
			xl=rtf;
			fl=f;
		} else {
			del=xh-rtf;
			xh=rtf;
			fh=f;
		}
		//dx=xh-xl;
		xm = (xl + xh) / 2.0;
		if (fabs(del) < xacc || f == 0.0) return rtf;
	}
	nrerror("Maximum number of iterations exceeded in rtmuller");
	return 0.0;
}
#undef MAXIT
