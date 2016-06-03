//
const double rho0 = 2300.0;
const double p0 = 1.0e5;
const double Gam_lq = 7.15;
const double c0 = 1500.0;
const double g = 10.0;
const double Kh = 4.33e-6;
const double kB = 1.38e-23;
const double Rg = 8.314;
const double T = 1000;
const double E0 = 5.1e-19;
const double mu0 =  0.003162278;
const double M_m = 3.0e-26; // Massa molecul
const double MH2O = 1.8e-2;
const double De = 2.3e-11; // 2.3e-11;
const double Pi = 3.1415;
const double sigma = 0.072;
const double Na = 6.0e23; // Число Авогадро

const double Rb_min = 1.0e-9; // м

//
double Xm, Zm, Z_max, Z_ch;
int Im,J_max, Ig;
int *Jm;
double hx, hz;
double Tm, tau, T_all, T_out, T_surf;
double E, J, J_e, kapp, V_m, d_m, N_0, V_0;
double rho_surf, Cp_temp, u_av, dMg, dP, P_eq, W_cr;

double *x, *z;
double **p;
double **Pg;
double **rho, **rho_n;
double **u, **u_n;
double **v, **v_n;
double **mu;
double **Cg;
double **Cp, **Cp0;
double **Nb;
double **Rb;
double **dRb, **ddRb;
double **Eps;
double **Mg;
double **D_eff;
double *GrX, *GrZ;
double **Nt;
int **Ind;


int i, j, l;
int out_num;
int i_out, j_out;

char out_name[25];
FILE *out_file, *cut_file;

int var_case;
int var_wall;
int var_Eq_state;
