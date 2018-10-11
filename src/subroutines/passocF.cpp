#include "GF_calls.h"
#include "Faddeeva.hh"
#include "rand_gsl.h"

double passocF(double r0, double tcurr, double Dtot, double bindrad, double alpha, double cof) {
	double fDt = 4.0 * Dtot * tcurr;
	double sq_fDt = sqrt(fDt);

	double f1 = cof * bindrad / r0;

	int i, j;
	double sep, dist;
	double sqrt_t = sqrt(tcurr);
	double a2 = alpha * alpha;
	double r1, term1, term2, e1, ef1, sum;
	double onemsirr;
	sep = (r0 - bindrad) / sq_fDt; //a

	e1 = 2.0 * sep * sqrt_t * alpha + a2 * tcurr;
	ef1 = sep + alpha * sqrt_t; //a+b
	term1 = erfc(sep);
	double ep1 = exp(e1);
	if (isinf(ep1)) {
		std::complex<double> z;
		z.real(0.0);// = 0.0;
		z.imag(ef1);// = ef1;
		
		//cout <<"Complex number: "<<z<<endl;
		std::complex<double> value;
		double relerr = 0;
		value = Faddeeva::w(z, relerr);
		double ea2 = exp(-sep * sep);
		term2 = ea2 * real(value);

	} else {
		term2 = exp(e1) * erfc(ef1);
	}
	sum = term1 - term2;
	sum *= f1;
	onemsirr = sum;
	return onemsirr; //1-sirr=passoc
	//  cout <<"s_irr: "<<sirr<<" time: "<<tcurr<<" unscaled: "<<term1-term2<<endl;
}
