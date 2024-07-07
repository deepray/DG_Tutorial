clc
clear all
close all


model     = 'Euler';
gas_const = 1.0;
gas_gamma = 1.4;
test_name = 'Sod';
rho_IC    =@(x) 1*(x<0.0) + 0.125*(x>=0.0);
vel_IC    =@(x) 0*x;
pre_IC    =@(x) 1*(x<0.0) + 0.1*(x>=0.0);

bnd_l     = -5;  
bnd_r     = 5;
mesh_pert = 0.0;
bc_cond   = {'D',1,'D',0.125;
             'D',0,'D',0.0;
             'D',1/(0.4),'D',0.1/(0.4)};  % For conserved variables
FinalTime = 2;
CFL       = 0.4;
K         = 200;
N         = 1;


Indicator      = 'NN';
ind_var        = 'prim';
Limiter        = 'NONE';
lim_var        = "char_stencil";
Visc_model     = "MDH";
c_A            = 2.0;
c_k            = 0.4;
c_max          = 0.5;

% Plot and save parameters
plot_iter  = 10;
ind_iter   = 2;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
var_ran    = [0,1.2; -0.2,1.5; 0,1.2];

% Call code driver
EulerDriver1D; 




