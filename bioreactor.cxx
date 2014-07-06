// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct Par {
	double ks;  
	double kd;  
	double eta;  
	double lambda_fi;  
	double mu_max;   
	double sr;   
	double dr;} p;

void Bioreactor( double v[], double dv[] ){
	double S, //substrate
	    	 X; //biomass
	     
	//Variables
	S = v[0];
	X = v[1];

	//The Model
	dv[0] = p.mu_max * ( S * X ) / ( p.ks + S ) - p.kd * X;
	dv[1] = p.sr * p.dr - p.dr * S - p.mu_max / p.eta * ( S * X ) / ( p.ks + S );
}

/* Runge Kutta integrator from numerical recipies plus improvements */
void rk4( void deri(double [], double []), double h[], int n, double dt ){
	#define naux 10

	int i;
	double k1[naux],k2[naux],k3[naux],k4[naux],h0[naux];
	double dt2, dt6;

	dt2=dt/2.;
	dt6=dt/6.;

	for (i = 0 ; i<n; i++)
		h0[i] = h[i];

	deri(h0,k1);
	for (i =0 ; i<n ; i++)
		h0[i]=h[i]+dt2*k1[i];

	deri(h0,k2);
	for (i =0 ; i<n ; i++)
		h0[i]=h[i]+dt2*k2[i];

	deri(h0,k3);
	for (i =0 ; i<n ; i++)
		h0[i]=h[i]+dt*k3[i];

	deri(h0,k4);
	for (i = 0 ; i<n ; i++)
		{h0[i]=h[i]+dt*k4[i];};

	for (i =0; i<n ; i++)
		h[i]=h[i]+dt6*(2.*(k2[i]+k3[i])+k1[i]+k4[i]);
}

//Main Function
int main( void )
{
	double v[2],
	t = 0,
	tau = 1000,
	dt = 0.1;
	v[0] = 0.0;
	v[1] = 100.0;
	p.ks = 1.0;  
	p.kd = 0.003;  
	p.eta = 1/0.32;  
  p.mu_max = 1.0;   
	p.sr = 100;   
	p.dr = 0.01;

	while (t < tau ){
		rk4( Bioreactor, v, 2, dt );
		t+=dt;
		printf("%4.2f %4.9f %4.9f\n", t, v[0], v[1]);
	}

	return 0;
}
