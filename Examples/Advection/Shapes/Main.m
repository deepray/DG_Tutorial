clc
clear all
close all


model     = 'Advection';
test_name = 'Shapes';
u_IC =@(x) 10*(x-0.2).*(x>=0.2).*(x<0.3)... 
            + 10*(0.4-x).*(x>=.3).*(x<0.4)...
            + 1*(x>=.6).*(x<0.8)...
            + 100*(x-1.0).*(1.2-x).*(x>=1.0).*(x<1.2); 

        
bnd_l     = 0.0;  
bnd_r     = 1.4;
mesh_pert = 0.0;
bc_cond   = {'P',0.0,'P',0.0};
FinalTime = 1.4;
CFL       = 0.2;
K         = 200;
N         = 4;


Indicator  = 'TVB'; TVBM = 10;
Limiter    = 'MINMOD';
Visc_model = 'NONE';


plot_iter  = 50;
ind_iter   = 50;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
var_ran    = [-0.2,1.5];

% Call code driver
ScalarDriver1D; 




