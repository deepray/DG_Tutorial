clc
clear all
close all


model     = 'Burgers';
test_name = 'Sine';
u_IC =@(x)  sin(2*pi*x);           

        
bnd_l     = 0.0;  
bnd_r     = 1.0;
mesh_pert = 0.0;
bc_cond   = {'P',0.0,'P',0.0};
FinalTime = 0.5;
CFL       = 0.2;
K         = 100;
N         = 4;


Indicator  = 'TVB'; TVBM = 10;
Limiter    = 'MINMOD';
Visc_model = 'NONE';


plot_iter  = 10;
ind_iter   = 5;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
var_ran    = [-1,1];

% Call code driver
ScalarDriver1D; 




