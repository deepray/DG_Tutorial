clc
clear all
close all


model     = 'Burgers';
test_name = 'Complex';
u_IC =@(x)  sin(1*pi*x).*(abs(x)>=1)...
            + 3*(x>=-1).*(x<-0.5)...
            + 1*(x>=-0.5).*(x<0.0)...
            + 3*(x>=0.0).*(x<0.5)...
            + 2*(x>=0.5).*(x<1.0);


bnd_l     = -4.0;  
bnd_r     = 4.0;
mesh_pert = 0.0;
bc_cond   = {'P',0.0,'P',0.0};
FinalTime = 0.4;
CFL       = 0.4;
K         = 200;
N         = 4;


Indicator  = 'MINMOD';
Limiter    = 'NONE';
Visc_model = 'MDH';
c_A        = 2.0;
c_k        = 0.4;
c_max      = 0.5;

plot_iter  = 40;
ind_iter   = 10;
ref_avail  = true;
ref_fname  = 'ref_soln.dat';
var_ran    = [-1.2,4];

% Call code driver
ScalarDriver1D; 




