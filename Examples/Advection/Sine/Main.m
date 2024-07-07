%CleanUp1D;
clc
clear all
close all


model     = 'Advection';
test_name = 'Sine';
u_IC =@(x) sin(10*pi*x); 


bnd_l     = 0;  
bnd_r     = 1.0;
mesh_pert = 0.0;
bc_cond   = {'P',0.0,'P',0.0};
FinalTime = 2;
CFL       = 0.6;
K         = 100;
N         = 4;


Indicator  = 'TVB'; TVBM = 10;
Limiter    = 'MINMOD';
Visc_model = 'NONE';


plot_iter  = 50;
ind_iter   = 20;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
var_ran    = [-1.2,1.5];

% Call code driver
ScalarDriver1D; 




