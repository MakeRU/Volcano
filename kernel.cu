#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <sstream>
#include <string>


#include "Data.h"


__device__ __host__ double StateEq(double rho, double Cg, double rho0, double c0, double Gam_lq, double p0, int var_Eq_state)
{
	double tmp;
		// return  c0*c0*(rho - rho0);
	switch (var_Eq_state)
	{
	case 0: { tmp = rho0*c0*c0 / Gam_lq * (pow(rho / (rho0*(1 - Cg)), Gam_lq) - 1.0) + p0;
		if (tmp < p0) tmp = p0; }; break;
	case 1: { tmp = rho0*c0*c0 * (rho / (rho0*(1 - Cg)) - 1.0) + p0;
		if (tmp < p0) tmp = p0; }; break;
	}

	return tmp;

	}

__device__ __host__ double StateEq1(double p, double Cg, double rho0, double c0, double Gam_lq, double p0, int var_Eq_state)
	{
		//return  rho0 + p / (c0*c0);
		switch (var_Eq_state)
		{
		case 0: {return rho0*(1 - Cg)*pow((p - p0)*Gam_lq / (rho0*c0*c0) + 1.0, 1.0 / Gam_lq); }; break;
		case 1: {return rho0*(1 - Cg)*( (p - p0)/ (rho0*c0*c0) + 1.0); }; break;
		}
	}

double FreeSurf(double T, double p)
{
//	double p_tmp;
//	if (T > T_surf) {return p0;}
//		else {return T* (p0-p) / (T_surf) + p;};
//	p_tmp = p / 1.2;
//	if (p_tmp < p0) {p_tmp = p0;};
	return p0;
}

int main()
{
	var_wall = 1;
	var_Eq_state = 1;
	out_num = 0;
	Tm = 0.0;
	T_all = 1.5;
	tau = 1.0e-6;
	T_out = 1000.0*tau;
	T_surf = 10.0*tau;
	hx = 0.25;
	hz = 1.0;
	Xm = 10.0;
	Zm = 1000.0;
	Z_max = 1.5*Zm;
	Z_ch = 7.5e3;

	Im = int(Xm / hx) + 1;
	J_max = int(Z_max / hz) + 1; 
	Jm = new int [Im];
	for (i = 0; i < Im; i++){
			Jm[i] = int(Zm/hz);
		}

	x = new double [Im];
	for (i = 0; i < Im; i++){
		x[i] = i*hx;
	}
	z = new double [J_max];
	for (j = 0; j < J_max; j++){
		z[j] = j*hz;
	}

	Ig = 5*Im;
	GrX = new double [Ig];
	GrZ = new double [Ig];
	for (i = 0; i < Ig; i++){
			GrX[i] = i * double (Xm / (Ig-1));
			GrZ[i] = Zm;
		}


	p = new double *[Im];
	Pg = new double *[Im];
	rho = new double *[Im];
	rho_n = new double *[Im];
	u = new double *[Im];
	v = new double *[Im];
	u_n = new double *[Im];
	v_n = new double *[Im];
	mu = new double *[Im];
	Cg = new double *[Im];
	Cp = new double *[Im];
	Cp0 = new double *[Im];
	Nb = new double *[Im];
	Rb = new double *[Im];
	dRb = new double *[Im];
	ddRb = new double *[Im];
	Eps = new double *[Im];
	Mg = new double *[Im];
	D_eff = new double *[Im];
	Nt = new double *[Im];
	Ind = new int *[Im];
	for (i = 0; i < Im; i++){
		p[i] = new double [J_max];
		Pg[i] = new double[J_max];
		rho[i] = new double [J_max];
		rho_n[i] = new double [J_max];
		u[i] = new double [J_max];
		v[i] = new double [J_max];
		u_n[i] = new double [J_max];
		v_n[i] = new double [J_max];
		mu[i] = new double [J_max];
		Cg[i] = new double [J_max];
		Cp[i] = new double [J_max];
		Cp0[i] = new double [J_max];
		Nb[i] = new double [J_max];
		Rb[i] = new double [J_max];
		dRb[i] = new double[J_max];
		ddRb[i] = new double[J_max];
		Eps[i] = new double [J_max];
		Mg[i] = new double [J_max];
		D_eff[i] = new double [J_max];
		Nt[i] = new double [J_max];
		Ind[i] = new int [J_max];
	}

	for (i = 0; i < Im; i++){
		for (j = 0; j < J_max; j++){
			p[i][j] = rho0 * g * (Z_ch - z[j]);	
			Pg[i][j] = p[i][j];
			Cp0[i][j] = Kh*sqrt(p[i][j]);
			Nb[i][j] = 0.0;// 1.0e13;
			Rb[i][j] = 0.0; // Rb_min;
			dRb[i][j] = 0.0;
			ddRb[i][j] = 0.0;
			Mg[i][j] = 0.0; // 4.0*Pi * Pg[i][j] *Rb[i][j]*Rb[i][j]*Rb[i][j]*MH2O / (3.0* Rg * T);
			Cp[i][j] = Cp0[i][j] - Nb[i][j] * Mg[i][j] / rho0;
			if (Cp[i][j] < 0) {Cp[i][j] = 1.0e-3;};
        	Cg[i][j] = Nb[i][j]* 4.0*Pi/3.0 *Rb[i][j]*Rb[i][j]*Rb[i][j];
        	//Cg[i][j] = Cg[i][j] / (1.0 + Cg[i][j]);
			rho[i][j] = StateEq1(p[i][j], Cg[i][j], rho0, c0, Gam_lq, p0, var_Eq_state);
       		E = E0*(1-12*Cp[i][j]);
       		mu[i][j] =mu0*exp(E/(kB*T));
		//	mu[i][j] = 1.0e7;
			u[i][j] = 0.0;
			v[i][j] = 0.0;
			u_n[i][j] = 0.0;
			v_n[i][j] = 0.0;
			Eps[i][j] = 0.0;
			D_eff[i][j] = 0.0;
		//	mu[i][j] = 1.0e7;
			Nt[i][j] = 0.0;
			Ind[i][j] = 0;

		}
	}
       
	rho_surf = rho[0][0];

	for (i = 0; i < Im; i++){
           	for (j = 1; j < J_max-1; j++){
           		if (j > Jm[i]) {
           			u[i][j] = 0.0;
           			v[i][j] = 0.0;
					mu[i][j] = 1.0;					
					p[i][j] = p0;
					rho[i][j] = StateEq1(p[i][j], Cg[i][j], rho0, c0, Gam_lq, p0, var_Eq_state);

              	}
            }
        }

	sprintf(out_name, "DataI/P");
	cut_file = fopen( out_name, "wt" );
	i_out = 0;
	j_out = 900;

	do
	{
		Tm = Tm + tau;
        printf("Time %5.6lf s \n", Tm);


        // ����� u

        for (i = 1; i < Im-1; i++){
        	for (j = 1; j < J_max-1; j++){
        		u_n[i][j] = u[i][j] + tau * ( 	-u[i][j] * (u[i+1][j] - u[i-1][j])/(2.0*hx)
        										-v[i][j] * (u[i][j+1] - u[i][j-1])/(2.0*hz) -
        										1.0/rho[i][j] * (p[i+1][j] - p[i][j])/hx +
        		 mu[i][j]/rho[i][j] * ( (u[i+1][j] -2.0 * u[i][j] +u[i-1][j])/(hx*hx) +
        				(u[i][j+1] -2.0 * u[i][j] +u[i][j-1])/(hz*hz))	+
						1 / rho[i][j] * ( (mu[i+1][j] - mu[i][j]) / hx * (u[i+1][j] - u[i-1][j]) / (2.0*hx) + 
						(mu[i][j + 1] - mu[i][j]) / hz + (u[i][j + 1] - u[i][j - 1]) / (2.0*hz)));
        	}
        }

        for (i = 1; i < Im-1; i++){
        	for (j = 1; j < J_max-1; j++){
        		u[i][j] = u_n[i][j];
        	}
        }
        // ��� ���������
       	for (j = 1; j < J_max-1; j++){
             u[0][j] = 0.0;
        }


        for (i = 1; i < Im-1; i++){
        	for (j = 1; j < J_max-1; j++){
				v_n[i][j] = v[i][j] + tau * (-10.0 - u[i][j] * (v[i + 1][j] - v[i - 1][j]) / (2.0*hx)
        										-v[i][j] * (v[i][j+1] - v[i][j-1])/(2.0*hz) -
        										1.0/rho[i][j] * (p[i][j+1] - p[i][j])/hz
        		+ mu[i][j]/rho[i][j] * ( (v[i+1][j] -2.0 * v[i][j] +v[i-1][j])/(hx*hx) +
        				(v[i][j+1] - 2.0 * v[i][j] + v[i][j-1])/(hz*hz))	+
						1 / rho[i][j] * ((mu[i + 1][j] - mu[i][j]) / hx * (v[i + 1][j] - v[i - 1][j]) / (2.0*hx) +
						(mu[i][j + 1] - mu[i][j]) / hz + (v[i][j + 1] - v[i][j - 1]) / (2.0*hz)));
        	}
        }

        for (i = 1; i < Im-1; i++){
        	for (j = 1; j < J_max-1; j++){
        		v[i][j] = v_n[i][j];
        	}
        }

        // ��� ���������
       	for (j = 1; j < J_max-1; j++){
             v[0][j] = v_n[1][j];
        }
		
		for (i = 0; i < Im; i++){
			v[i][0] = v[i][1];
		}

       	if (var_wall == 1) {
          for (j = 1; j < J_max-1; j++){
                 u_av = 0.0;
            for (i = 0; i < Im-1; i++){
            	 u_av =  u_av + v[i][j];
            }
            u_av = u_av / (Im-1);
			if (Cg[Im-2][j] < 0.4 ) {v[Im-1][j] = u_av * Cg[Im-2][j] / (1.0 - Cg[Im-2][j]);}
			else  {v[Im-1][j] = u_av * Cg[Im-2][j];};
  	      }
       	}

        // ����� rho

        for (i = 1; i < Im-1; i++){
          	for (j = 1; j < J_max-1; j++){
          		rho_n[i][j] = rho[i][j] + tau * (-u[i][j] * (rho[i+1][j] - rho[i-1][j])/(2.0*hx)
          										 -v[i][j] * (rho[i][j+1] - rho[i][j-1])/(2.0*hz) -
          					rho[i][j] * ((u[i][j] - u[i-1][j])/hx + (v[i][j] - v[i][j-1])/hz ));
				if (rho_n[i][j] < rho0 *0.25) { rho_n[i][j] = rho0 *0.25; };
          	}
          }

        // ��� ���������  � �������
       	for (j = 1; j < J_max-1; j++){
             rho_n[0][j] = rho_n[1][j];
             rho_n[Im-1][j] = rho_n[Im-2][j];
        }

        // ����� p

        for (i = 1; i < Im-1; i++){
        	for (j = 1; j < J_max-1; j++){
        		rho[i][j] = rho_n[i][j];
				p[i][j] = StateEq(rho[i][j], Cg[i][j], rho0, c0, Gam_lq, p0, var_Eq_state);
        	}
        }

		        // ��� ���������  � �������
       	for (j = 1; j < J_max-1; j++){
             p[0][j] = p[1][j];
			 rho[0][j] = StateEq1(p[0][j], Cg[0][j], rho0, c0, Gam_lq, p0, var_Eq_state);
             p[Im-1][j] = p[Im-2][j];
			 rho[Im - 1][j] = StateEq1(p[Im - 1][j], Cg[Im - 1][j], rho0, c0, Gam_lq, p0, var_Eq_state);
        }

		// Nucleation

		for (i = 0; i < Im; i++){
			for (j = 1; j <= Jm[i]; j++){
				switch(Ind[i][j])
				{ 
					// Nucleation
				case 0 : {		
					P_eq = Cp[i][j] * Cp[i][j] / (Kh*Kh);
					dP = P_eq - p[i][j];
					if (dP > 0.0) 
					{
						W_cr = 16.0*Pi*sigma*sigma*sigma / (3.0*dP*dP);
						kapp = 3.0;
						V_0 = 1.0e-6;
						N_0 = rho0 * Cp[i][j] / M_m;
						V_m = 0.018 / (1000 * Na);
						d_m = pow(6.0 / (Pi * N_0), 1.0 / 3.0);
						J_e = 2.0*N_0*N_0*V_m*De / d_m * sqrt(sigma / (kB*T));
						J = J_e*exp(-W_cr / (kB*T));
						Nt[i][j] = Nt[i][j] + V_0*J*tau;
						if (Nt[i][j] > 1.0) 
						{
							Nb[i][j] = 0.6* pow((kapp*kapp*kapp - 1.0) / (kapp*kapp*kapp*kapp*kapp) * J / De, 3.0 / 5.0);
							Rb[i][j] = 2.0*sigma / (dP);
							Pg[i][j] = p[i][j];
							Mg[i][j] = 4.0*Pi * Pg[i][j] * Rb[i][j] * Rb[i][j] * Rb[i][j] * MH2O / (3.0* Rg * T);
							Cp[i][j] = Cp0[i][j] - Nb[i][j] * Mg[i][j] / rho0;
							Cg[i][j] = Nb[i][j] * 4.0*Pi / 3.0 *Rb[i][j] * Rb[i][j] * Rb[i][j];
							//Cg[i][j] = Cg[i][j] / (1.0 + Cg[i][j]);
							E = E0*(1 - 12 * Cp[i][j]);
							mu[i][j] = mu0*exp(E / (kB*T));
						//	mu[i][j] = 1.0e7;
							rho[i][j] = StateEq1(p[i][j], Cg[i][j], rho0, c0, Gam_lq, p0, var_Eq_state);
							Ind[i][j] = 1;
						}
				}
				}
				break;
					// Grows bubbles 
				case 1: {
					dRb[i][j] = Rb[i][j] * (Pg[i][j] - p[i][j]) / (4.0* mu[i][j]);
					Rb[i][j] = Rb[i][j] + dRb[i][j] * tau;
					if (Rb[i][j] < Rb_min) { Rb[i][j] = Rb_min; };
					if (p[i][j] > p0) { Cp_temp = Kh*sqrt(p[i][j]); }
					else { Cp_temp = Kh*sqrt(p0); };
					dMg = 4.0*Pi* Rb[i][j] * rho0 * De * (Cp[i][j] - Cp_temp) *tau;
					if ((dMg > 0.0) && (dMg*Nb[i][j] > Cp[i][j] * rho0)) { dMg = Cp[i][j] * rho0 / Nb[i][j]; };
					//	if ((dMg < 0.0) && (dMg > Cp[i][j] * rho0)) { dMg = Cp[i][j] * rho0 / Nb[i][j]; };
					Mg[i][j] = Mg[i][j] + dMg;
					if (Mg[i][j] < 0.0) { Mg[i][j] = 4.0*Pi * p[i][j] * Rb[i][j] * Rb[i][j] * Rb[i][j] * MH2O / (3.0* Rg * T); };
					Pg[i][j] = 3.0* Mg[i][j] * Rg * T / (4.0*Pi *Rb[i][j] * Rb[i][j] * Rb[i][j] * MH2O);
					Cp[i][j] = Cp0[i][j] - Nb[i][j] * Mg[i][j] / rho0;
					if (Cp[i][j] < 0) { Cp[i][j] = Kh*sqrt(p0); };
					E = E0*(1 - 12 * Cp[i][j]);
					mu[i][j] = mu0*exp(E / (kB*T));			
				//	mu[i][j] = 1.0e7;
					Cg[i][j] = Nb[i][j] * 4.0*Pi / 3.0 *Rb[i][j] * Rb[i][j] * Rb[i][j];
					Cg[i][j] = Cg[i][j] / (1.0 + Cg[i][j]);
				//	if (Cg[i][j] > 0.75) { Cg[i][j] = 0.75; };
					p[i][j] = StateEq(rho[i][j], Cg[i][j], rho0, c0, Gam_lq, p0, var_Eq_state);
				}
						break;
				case 2: {
					if (p[i][j] > p0) { Eps[i][j] = rho0 / (MH2O * p[i][j] / (Rg * T))* (Cp[i][j] - Kh*sqrt(p[i][j])); }
					else { Eps[i][j] = rho0 / (MH2O * p0 / (Rg * T))* (Cp[i][j] - Kh*sqrt(p0)); };
					if (Eps[i][j] < 0.0) { Eps[i][j] = 0.0; };
					if (Eps[i][j] > 2.0) { D_eff[i][j] = De*12.0 / Pi * Eps[i][j] * Eps[i][j]; }
					else { D_eff[i][j] = De * 2.0 * Eps[i][j]; };
					Rb[i][j] = Rb[i][j] + tau*D_eff[i][j] / (2.0*Rb[i][j]);
					if (Rb[i][j] < Rb_min) { Rb[i][j] = Rb_min; };
					if (p[i][j] > p0) { Pg[i][j] = p[i][j]; }
					else {
						Pg[i][j] = p0;
					}
					Mg[i][j] = (MH2O * Pg[i][j] / (Rg * T)) * 4.0*Pi / 3.0 *Rb[i][j] * Rb[i][j] * Rb[i][j]; 
					Cp[i][j] = Cp0[i][j] - Nb[i][j] * Mg[i][j] / rho0;
					if (Cp[i][j] < 0) { Cp[i][j] = Kh*sqrt(p0); };
					Cg[i][j] = Nb[i][j] * 4.0*Pi / 3.0 *Rb[i][j] * Rb[i][j] * Rb[i][j];
					Cg[i][j] = Cg[i][j] / (1.0 + Cg[i][j]);
					p[i][j] = StateEq(rho[i][j], Cg[i][j], rho0, c0, Gam_lq, p0, var_Eq_state);
					E = E0*(1 - 12 * Cp[i][j]);
					mu[i][j] = mu0*exp(E / (kB*T));
				} 
						break;
				case 3: {
					ddRb[i][j] = (Pg[i][j] - p[i][j]) / (Rb[i][j] * rho[i][j]) -
						(4.0* mu[i][j] * dRb[i][j]) / (Rb[i][j] * rho[i][j] * Rb[i][j]) -
						3.0*dRb[i][j] * dRb[i][j] / (2.0*Rb[i][j]);
					dRb[i][j] = dRb[i][j] + ddRb[i][j] * tau;
					Rb[i][j] = Rb[i][j] + dRb[i][j] * tau;
					if (Rb[i][j] < Rb_min) { Rb[i][j] = Rb_min; };
					if (p[i][j] > p0) { Cp_temp = Kh*sqrt(p[i][j]); }
					else { Cp_temp = Kh*sqrt(p0); };
					dMg = 4.0*Pi* Rb[i][j] * rho0 * De * (Cp[i][j] - Cp_temp) *tau;
					if ((dMg > 0.0) && (dMg*Nb[i][j] > Cp[i][j] * rho0)) { dMg = Cp[i][j] * rho0 / Nb[i][j]; };
					//	if ((dMg < 0.0) && (dMg > Cp[i][j] * rho0)) { dMg = Cp[i][j] * rho0 / Nb[i][j]; };
					Mg[i][j] = Mg[i][j] + dMg;
					if (Mg[i][j] < 0.0) { Mg[i][j] = 4.0*Pi * p[i][j] * Rb[i][j] * Rb[i][j] * Rb[i][j] * MH2O / (3.0* Rg * T); };
					Pg[i][j] = 3.0* Mg[i][j] * Rg * T / (4.0*Pi *Rb[i][j] * Rb[i][j] * Rb[i][j] * MH2O);
					Cp[i][j] = Cp0[i][j] - Nb[i][j] * Mg[i][j] / rho0;
					if (Cp[i][j] < 0) { Cp[i][j] = Kh*sqrt(p0); };
					E = E0*(1 - 12 * Cp[i][j]);
					mu[i][j] = mu0*exp(E / (kB*T));
					Cg[i][j] = Nb[i][j] * 4.0*Pi / 3.0 *Rb[i][j] * Rb[i][j] * Rb[i][j];
					Cg[i][j] = Cg[i][j] / (1.0 + Cg[i][j]);
					//	if (Cg[i][j] > 0.75) { Cg[i][j] = 0.75; };
					p[i][j] = StateEq(rho[i][j], Cg[i][j], rho0, c0, Gam_lq, p0, var_Eq_state);
				} break;

				}
		
			}
		}


 /*     	for (i = 0; i < Im; i++){
			for (j = 1; j <= Jm[i]; j++){
				dRb[i][j] = Rb[i][j] * (Pg[i][j] - p[i][j]) / (4.0* mu[i][j]);
				Rb[i][j] = Rb[i][j] + dRb[i][j] * tau;
				if (Rb[i][j] < Rb_min) { Rb[i][j] = Rb_min; };
				if (p[i][j] > p0) {Cp_temp = Kh*sqrt(p[i][j]);}
				else {Cp_temp = Kh*sqrt(p0);};
				dMg = 4.0*Pi* Rb[i][j] * rho0 * De * (Cp[i][j] - Cp_temp) *tau;
				if ((dMg > 0.0) && (dMg*Nb[i][j] > Cp[i][j] * rho0)) { dMg = Cp[i][j] * rho0 / Nb[i][j]; };
			//	if ((dMg < 0.0) && (dMg > Cp[i][j] * rho0)) { dMg = Cp[i][j] * rho0 / Nb[i][j]; };
				Mg[i][j]= Mg[i][j] + dMg;
				if (Mg[i][j] < 0.0) { Mg[i][j] = 4.0*Pi * p[i][j] * Rb[i][j] * Rb[i][j] * Rb[i][j] * MH2O / (3.0* Rg * T); };
				Pg[i][j] = 3.0* Mg[i][j] * Rg * T / (4.0*Pi *Rb[i][j] * Rb[i][j] * Rb[i][j]* MH2O);
				Cp[i][j] = Cp0[i][j] - Nb[i][j]* Mg[i][j]/rho0;
				if (Cp[i][j] < 0) {Cp[i][j] = Kh*sqrt(p0);};
				E = E0*(1-12*Cp[i][j]);
				mu[i][j] =mu0*exp(E/(kB*T));
				Cg[i][j] = Nb[i][j]* 4.0*Pi/3.0 *Rb[i][j]*Rb[i][j]*Rb[i][j];
				Cg[i][j] = Cg[i][j] / (1.0 + Cg[i][j]);
				if (Cg[i][j] > 0.75) { Cg[i][j] = 0.75; };
				p[i][j] = StateEq(rho[i][j], Cg[i][j], rho0, c0, Gam_lq, p0);

			}
		}
		*/

  /*    //     	for (i = 0; i < Im; i++){
     	//		for (j = 1; j <= Jm[i]; j++){
       	{{ i=i_out; j = j_out;
     				ddRb[i][j] = (Pg[i][j] - p[i][j]) / (Rb[i][j] * rho[i][j]) -
     						(4.0* mu[i][j]*dRb[i][j])/(Rb[i][j] * rho[i][j]*Rb[i][j]) -
     								3.0*dRb[i][j]*dRb[i][j]/(2.0*Rb[i][j]);
     				dRb[i][j] = dRb[i][j] + ddRb[i][j] * tau;
     				Rb[i][j] = Rb[i][j] + dRb[i][j] * tau;
     				if (Rb[i][j] < Rb_min) { Rb[i][j] = Rb_min; };
     				Pg[i][j] = Mg[i][j] / (4.0*Pi / 3.0 *Rb[i][j] * Rb[i][j] * Rb[i][j]) * Rg * T / MH2O;
     				Cp[i][j] = Cp0[i][j] - Nb[i][j]* Mg[i][j]/rho0;
     				if (Cp[i][j] < 0) {Cp[i][j] = Kh*sqrt(p0);};
     				Cg[i][j] = Nb[i][j]* 4.0*Pi/3.0 *Rb[i][j]*Rb[i][j]*Rb[i][j];
     				Cg[i][j] = Cg[i][j] / (1.0 + Cg[i][j]);
     				p[i][j] = StateEq(rho[i][j], Cg[i][j], rho0, c0, Gam_lq, p0);

     			}
     		}
*/
 /*       for (i = 0; i < Im; i++){
        	for (j = 1; j <= Jm[i]; j++){
        		if (p[i][j] > p0) {Eps[i][j] = rho0/ (MH2O * p[i][j] / (Rg * T))* (Cp[i][j] - Kh*sqrt(p[i][j]));}
        		else {Eps[i][j] = rho0/ (MH2O * p0 / (Rg * T))* (Cp[i][j] - Kh*sqrt(p0));};
				if (Eps[i][j] < 0.0 ) {Eps[i][j] = 0.0;};
           		if (Eps[i][j] > 2.0 ) {D_eff[i][j] = De*12.0/Pi * Eps[i][j] * Eps[i][j];}
           		else {D_eff[i][j] = De * 2.0 * Eps[i][j];};
        		Rb[i][j] = Rb[i][j] + tau*D_eff[i][j] / (2.0*Rb[i][j]);
        		if (Rb[i][j] < Rb_min) {Rb[i][j] = Rb_min;};
				if (p[i][j] > p0) {Mg[i][j] = (MH2O * p[i][j] / (Rg * T)) * 4.0*Pi/3.0 *Rb[i][j]*Rb[i][j]*Rb[i][j];}
				else {Mg[i][j] = (MH2O * p0 / (Rg * T)) * 4.0*Pi/3.0 *Rb[i][j]*Rb[i][j]*Rb[i][j];}
        		Cp[i][j] = Cp0[i][j] - Nb[i][j]* Mg[i][j]/rho0;
				if (Cp[i][j] < 0) {Cp[i][j] = Kh*sqrt(p0);};
        		Cg[i][j] = Nb[i][j]* 4.0*Pi/3.0 *Rb[i][j]*Rb[i][j]*Rb[i][j];
        		Cg[i][j] = Cg[i][j] / (1.0 + Cg[i][j]);
        		p[i][j] = StateEq(rho[i][j], Cg[i][j], rho0, c0, Gam_lq, p0);
        		E = E0*(1-12*Cp[i][j]);
        		mu[i][j] =mu0*exp(E/(kB*T));
        	}
        }
*/
            // ��� ��������� ������������
		
        GrZ[Ig-1] = GrZ[Ig-1] + v[Im-1][Jm[Im-1]]*tau;
        Jm[Im-1] = int (GrZ[Ig-1] / hz);
		for (l = Ig-2; l >= 0; l--){
			i= int (GrX[l] /hx);
        //  GrX[l] = GrX[l] + u[i][Jm[i]]*tau;
          GrZ[l] = GrZ[l] + v[i][Jm[i]]*tau;
		if (GrZ[l] > Z_max) {GrZ[l] = Z_max;};
		  if (GrZ[l] < GrZ[l+1]) {GrZ[l] = GrZ[l+1];};
          Jm[i] = int (GrZ[l] / hz);
		}


        for (i = 0; i < Im; i++){
           	for (j = 1; j < J_max; j++){
           		if (j > Jm[i]) {
           			u[i][j] = u[i][Jm[i]];
           			v[i][j] = v[i][Jm[i]];
					Nb[i][j] = Nb[i][Jm[i]];
					Rb[i][j] = Rb[i][Jm[i]];
           			Cp[i][j] = Cp[i][Jm[i]];
           			Cg[i][j] = Cg[i][Jm[i]];
					Mg[i][j] = Mg[i][Jm[i]];
					Pg[i][j] = Pg[i][Jm[i]];
					mu[i][j] = mu[i][Jm[i]];
					Ind[i][j] = Ind[i][Jm[i]];
					
			/*		u[i][j] = u[i][1000];
           			v[i][j] = v[i][1000];
           			Rb[i][j] = Rb[i][1000];
           			Cp[i][j] = Cp[i][1000];
           			Cg[i][j] = Cg[i][1000];
					mu[i][j] = mu[i][1000]; */
           			p[i][j] = FreeSurf(Tm, p[i][j]);
					rho[i][j] = StateEq1(p[i][j], Cg[i][j], rho0, c0, Gam_lq, p0, var_Eq_state);
              	}
            }
        }





        //

	//	if (Tm > 0.035) {T_out = tau; }
        // ������ ������
     //   if ( (Tm > out_num * T_out) || (Tm > 57.0e-3)){
		if (Tm > out_num * T_out) {

        if (out_num < 100000) {sprintf(out_name, "Data/%d", out_num);};
        if (out_num < 10000) {sprintf(out_name, "Data/0%d", out_num);};
        if (out_num < 1000) {sprintf(out_name, "Data/00%d", out_num);};
        if (out_num < 100) {sprintf(out_name, "Data/000%d", out_num);};
        if (out_num < 10) {sprintf(out_name, "Data/0000%d", out_num);};
        out_file = fopen( out_name, "wt" );
    	for (i = 0; i < Im; i++){
    		for (j = 0; j <= Jm[i]; j++){
    			fprintf( out_file, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \n",
    					x[i], z[j], p[i][j]/p0, rho[i][j], Cg[i][j], mu[i][j], u[i][j], v[i][j], Cp[i][j], Rb[i][j]*1.0e6, dRb[i][j], Pg[i][j]/p0, Nb[i][j]);
    		}
    		fprintf( out_file, "\n");
    	}
    	fclose( out_file );

    	if (out_num < 100000) {sprintf(out_name, "DataG/%d", out_num);};
    	if (out_num < 10000) {sprintf(out_name, "DataG/0%d", out_num);};
    	if (out_num < 1000) {sprintf(out_name, "DataG/00%d", out_num);};
    	if (out_num < 100) {sprintf(out_name, "DataG/000%d", out_num);};
    	if (out_num < 10) {sprintf(out_name, "DataG/0000%d", out_num);};
    	out_file = fopen( out_name, "wt" );
    	for (l = 0; l < Ig; l++){
			i= int (GrX[l] /hx);
    		fprintf( out_file, "%10.8lf \t %10.8lf \t %d \n",GrX[l], GrZ[l], Jm[i]);
    		}
    	fclose( out_file );

    	if (out_num < 100000) {sprintf(out_name, "DataB/%d", out_num);};
    	if (out_num < 10000) {sprintf(out_name, "DataB/0%d", out_num);};
    	if (out_num < 1000) {sprintf(out_name, "DataB/00%d", out_num);};
    	if (out_num < 100) {sprintf(out_name, "DataB/000%d", out_num);};
    	if (out_num < 10) {sprintf(out_name, "DataB/0000%d", out_num);};
    	out_file = fopen( out_name, "wt" );
    	for (j = 0; j <= Jm[0]; j++){
			fprintf(out_file, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \n", 
				z[j], p[0][j] / p0, v[0][j], Nb[0][j], Cp[0][j], Cg[0][j]);
    		}
    	fclose( out_file );


    	out_num = out_num + 1;
      	

    	fprintf( cut_file, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %e \t %e \t %e \t %e \t %lf \t %e\n",
    	    					Tm*1.0e3, p[i_out][j_out]/p0, mu[i_out][j_out], Cp[i_out][j_out],
    	    					Eps[i_out][j_out],D_eff[i_out][j_out]/De, Rb[i_out][j_out]*1.0e6, Mg[i_out][j_out],
    	    					Cg[i_out][j_out], dRb[i_out][j_out], Pg[i_out][j_out]/p0, Nt[i_out][j_out]);
	  }							
	}  while (Tm < T_all);

	fclose( cut_file );
}
