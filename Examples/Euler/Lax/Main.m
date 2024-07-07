clc
clear all
close all


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'Lax';
rho_IC =@(x) 0.445*(x<0) + 0.5*(x>=0.0);
vel_IC =@(x) 0.698*(x<0) + 0.0*(x>=0.0);
pre_IC =@(x) 3.528*(x<0) + 0.571*(x>=0.0);



bnd_l     = -5.0;  
bnd_r     = 5.0;
mesh_pert = 0.1;
bc_cond   = {'D',0.445,'D',0.5;
             'D',0.445*0.698,'N',0.0;
             'D',0.5*0.445*0.698^2 + 3.528/0.4,'D',0.571/0.4};
FinalTime = 1.3;
CFL       = 0.4;
K         = 200;
N         = 1;


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
var_ran    = [0,1.5; -1,2.5; 0,4];

% Call code driver
EulerDriver1D; 







