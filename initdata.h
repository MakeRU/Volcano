void InitVolcano()
{
	for (i = 0; i < Im; i++){
		for (j = 0; j < J_max; j++){
			Cg[i][j] = 0.0;
			p[i][j] = rho0 * g * (Z_ch - z[j]);
			rho[i][j] = StateEq1(p[i][j], Cg[i][j], rho0, c0, Gam_lq, p0);
			Cp0[i][j] = Kh*sqrt(p[i][j]);
			E = E0*(1-12*Cp0[i][j]);
			mu[i][j] =mu0*exp(E/(kB*T));
			u[i][j] = 0.0;
			v[i][j] = 0.0;
			u_n[i][j] = 0.0;
			v_n[i][j] = 0.0;
			Cp[i][j] = Cp0[i][j];
			Nb[i][j] = 1.0e12;
			Rb[i][j] = 1.0e-9;
			Eps[i][j] = 0.0;
			Mg[i][j] = 0.0;
			D_eff[i][j] = 0.0;
			mu[i][j] = 10.0;

		}
	}
 	
	for (i = 0; i < Im; i++){
			Jm[i] = int(Zm/hz);
		}

	for (i = 0; i < Im; i++){
           	for (j = 1; j < J_max-1; j++){
           		if (j > Jm[i]) {
           			u[i][j] = 0.0;
           			v[i][j] = 0.0;
					mu[i][j] = 10.0;
           			rho[i][j] = rho0;
           			p[i][j] = StateEq(rho[i][j], Cg[i][j], rho0, c0, Gam_lq, p0);
              	}
            }
        }


}

void InitBreak()
{
	//var_case = 5; // 5-Break X  6-Break Z
		for (i = 0; i < Im; i++){
			Jm[i] = int(Z_max/hz);
		}

	for (i = 0; i < Im; i++){
		for (j = 0; j < J_max; j++){
			Cg[i][j] = 0.0;
			switch (var_case) {
			case 5:	if (i > Im / 2) {p[i][j] = p0;} 
					else {p[i][j] = 1700.0*p0;}; 
					break;
			case 6:	if (j > J_max / 2) {p[i][j] = p0;} 
					else {p[i][j] = 1700.0*p0;}; 
					break;
			}
			rho[i][j] = StateEq1(p[i][j], Cg[i][j], rho0, c0, Gam_lq, p0);
			Cp0[i][j] = Kh*sqrt(p[i][j]);
			E = E0*(1-12*Cp0[i][j]);
			mu[i][j] = 0.0;
			u[i][j] = 0.0;
			v[i][j] = 0.0;
			u_n[i][j] = 0.0;
			v_n[i][j] = 0.0;
			Cp[i][j] = Cp0[i][j];
			Nb[i][j] = 1.0e12;
			Rb[i][j] = 1.0e-9;
			Eps[i][j] = 0.0;
			Mg[i][j] = 0.0;
			D_eff[i][j] = 0.0;
		//	mu[i][j] = 10.0;

		}
	}
       

}
