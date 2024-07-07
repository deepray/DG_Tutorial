clc
clear all
close all


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'LowDen';
rho_IC =@(x) 0*x + 7.0;
vel_IC =@(x) (x<0)*(-1.0)+ (x>=0.0)*1.0;
pre_IC =@(x) 0*x + 0.2;


bnd_l     = -1.0;  
bnd_r     = 1.0;
mesh_pert = 0.1;
bc_cond   = {'D',7,'D',7;
             'D',-7,'N',7.0;
             'D',4,'D',4.0};
FinalTime = 0.6;
CFL       = 0.4;
K     = 128;
N         = 4;


Indicator = 'NN';
ind_var   = 'prim';
Limiter   = 'MINMOD';
lim_var   = "char_stencil";
Visc_model     = "NONE";


% Plot and save parameters
plot_iter  = 10;
ind_iter   = 2;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
var_ran    = [0,9; -1.5,2.5; 0,0.3];

% Call code driver
EulerDriver1D; 









