# Discontinous Galerkin schemes for conservation laws (in 1D)

###Created by: Deep Ray, University of Maryland (deepray@umd.edu)
###Webpage: deepray.github.io 
###Date : 2 July, 2024

This repository contains tutoral MATLAB code for RKDG-solvers to solve 1D conservation laws. It is based on the code available in [DGANN](https://bitbucket.com/deepray/dgann) and [DGANN-AV](https://github.com/nickdisca/DGANN_AV). 

### Useful references to understand the code 
1. *Nodal Disontinuous Galerkin methods*; Hesthaven and Warburton. 
2. *An artificial neural network as a troubled-cell indicator*; Ray and Hesthaven. [[article]](https://www.sciencedirect.com/science/article/pii/S0021999118302547) [[preprint]](https://infoscience.epfl.ch/record/232425/files/manuscript_revised.pdf?version=1). 
3. *Controlling oscillations in high-order Discontinuous Galerkin schemes using artificial viscosity tuned by neural networks*; Discacciati, Hesthaven, and Ray [[article]](https://doi.org/10.1016/j.jcp.2020.109304) [[preprint]](https://infoscience.epfl.ch/record/263616?ln=en).
4. *Controlling oscillations in high-order schemes using neural networks*; Masters thesis by Discacciati [[thesis]](https://infoscience.epfl.ch/record/263615?ln=en&v=pdf).

##Table of contents 

* [Running the code](#markdown-header-running-the-code)
  * [Scalar 1D](#markdown-header-scalar-1d)
  * [Shallow water equations 1D](#markdown-header-shallow-water-1d)
  * [Euler equations 1D](#markdown-header-euler-1d)


##Running the code 
After cloning the git repository, execute **mypath.m** from the parent directory in MATLAB. This will set all the neccesary paths to use the solver. The various test cases need to be run from the **Examples** directory or its sub-directories. Currently, the 1D solver supports linear advection, Burgers' equation, Buckley-Leverett, and the compressible Euler equations.

####Scalar 1D
The basic structure of the example script is as follows.

~~~matlab
clc
clear all
close all


model = 'Advection';
test_name = 'Sine';
u_IC =@(x) sin(10*pi*x);   

bnd_l     = 0.0;  
bnd_r     = 1.4;
mesh_pert = 0.0;
bc_cond   = {'P',0.0,'P',0.0};
FinalTime = 1.4;
CFL       = 0.2;
K         = 100;
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
rk_comb    = false;
var_ran    = [-1.2,1.5];

% Call code driver
ScalarDriver1D; 
~~~

* The `model` flag sets the type of scalar model which is being solved. The following scalar models are currently available:
 * `'Advection'`: Linear advection equation with the advection speed set to 1.
 * `'Burgers'` : Burgers equation with the flux u^2/2.
 * `'BuckLev'`: Buckley-Leverett equation with flux contant  set to 0.5
* `test_name` is used to declare the name of the test. This is used for creating various save files.
* `u_IC` is used to set the initial condition for the problem.
* `bnd_l` and `bnd_r` define the left and right domain boundaries, which is descritized using `Nelem` number of elements. The degree of the polynomial in each element/cell is set using `N`. 
* `mesh_pert` is used to randomly perturb the interior cell-interfaces starting from a unifrom mesh.
* `bc_cond` is used to set the left and right boundary conditions. It is a cell of the form
 `{LEFT_BC_TYPE, LEFT_BC_VAL, RIGHT_BC_TYPE,RIGHT_BC_VAL}`. The `BC_TYPES` can be set as:
 * `'P'`: Periodic boundary conditions. In this case, both boundary type must be set to `'P'`.
 * `'N'`: Neumann boundary condition.
 * `'D'`: Dirichlet boundary condition, with the imposed Dirichlet value given by `BC_VAL`. Note that `BC_VAL` is ignored if `BC_TYPE` is not set to `'D'`. 
* The final simulation time is set using `FinalTime`, while the time step is chosen using `CFL`.
* `K` is the number of elements/cells in the mesh.
* `N` sets the order of the basis.
* The troubled-cell indicator is set using `Indicator`. The following options are currently available:
 * `'NONE'`: No cells are flagged.
 * `'ALL'`: All cells are flagged.
 * `'MINMOD'`: Uses the basic minmod limiter.
 * `'TVB'`: Uses the modified minmod-type TVB limiter. If this is chosen, then one also needs to set the variable `TVBM` to a positive number. Note that if this constant is set to 0, then the limiter reduces to the usual minmod limiter.
 * `'NN'`: Uses the trained neural network. 
* The limiter usd to reconstruct the solution in each troubled-cell,is set using `Limiter`. The following options are currently available:
 * `'NONE'`: No limiting is applied.
 * `'MINMOD'`: MUSCL reconstruction using minmod limiter.
* The artificial viscosity model is set using `Visc_model`. **NOTE** that you cannot use a viscosity model and limiting simultaneously. The following options are currently available:
 * `'NONE'`: No artificial viscosity is added.
 * `'EV'`: Uses the entropy viscosity model. If this is chosen, two nonnegative parameters `c_E` and `c_max` have to be specified.
 * `'MDH'`: Uses the highest modal decay model. If this is chosen, three nonnegative parameters `c_A`,`c_k`, and `c_max` have to be specified.
 * `'MDA'`: Uses the averaged modal decay model. If this is chosen, a nonnegative parameter `c_max` has to be specified.
* The flag `plot_iter` is used for visualization purposes. The solution plots are shown after every `plot_iter` number of iterations during the simulation. 
* `ind_iter` controls the frequency with which the time-history of the flagged troubled-cells or the spatial viscosity values are saved to file.
* If a reference/exact solution data file is available, then set `ref_avail` to `true` and the (relative) file name of the referene solution in `ref_fname`. 
* `var_ran` is used to set the ylim for the solution plot. This should be a array of size (1,2).
* The main driver script `ScalarDriver1D` is called once all the flags have been set.

[back to table of contents](#markdown-header-table-of-contents)
 
####Euler 1D 
The basic structure of the example script is as follows.

~~~matlab
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


Indicator      = 'MINMOD';
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
~~~

Most of the structure is similar to the Scalar 1D script. The differences are described below. 

* The `model` needs to be set as `'SWE'`. This should not be changed.
* The gas constant and ratio of specific heats is set using `gas_const` and `gas_gamma`.
* The initial density, velocity and pressure are set using `rho_IC`, `vel_IC` and `pre_IC`.
* The `bc_cond` has the same format as earlier, although now it has three rows of parameters. The first row corresponds to density, the second corresponds to the momentum, and the third to energy. Note that the boundary condition are set for the conserved variables.
* For systems of conservation laws, there are various choices for the variables to be used for troubled-cell detection. For the Euler equations, this choice is made via the flag `ind_var`, with the following options (the troubled-cells flagged for each variable is pooled together):
 * `'density'`: Only the density is used
 * `'velocity'`: Only the velocity is used
 * `'pressure'`: Only the pressure is used
 * `'prim'`: The primitive variables i.e., density, velocity and pressure, are used.
 * `'con'`: The conserved variables i.e., density, momentum and energy, are used. 
* As was the case with detection, there are several options for the variables which can be reconstructed. This is set using the flag `lim_var`, with the following options:
 * `'prim'`:  The primitive variables i.e., depth and velocity, are recontructed.
 * `'con'`:  The conserved variables i.e., depth and discharge, are recontructed.
 * `'char_cell'`: The local characterictic variables are reconstructed. These are obtained cell-by-cell using the linearized transformation operators. More precisely, the transformation matrix in each cell is evaluated using the cell-average value, following which the conserved variables are transformed to the characteristic variables in that cell. The same transformation is used to retrieve the conserved variables after limiting the characteristic variables. 
 * `'char_stencil'`: The local characterictic variables obtained cell-by-cell can introduce spurious oscillations in the solution. One can also obtain the local characteristic variables, stencil-by stencil. More precisely, for a given reconstruction stencil of 3-cells, the transformation matrix is evaluated using the cell-average value of the central cell, following which the conserved variables are transformed to the characteristic variables in every cell of that stencil. The transformed variables are used to obtain the reconstructed characteristic variables in the central cell. Note that this approach can be 3 times more expensive than the `'char_cell'` approach.
* `var_ran` is used to set the ylim for the solution plots, with the format `[rho_min,rho_max ; velocity_min, velocity_max ; pressure_min, pressure_max]`.
* The main driver script `EulerDriver1D` is called once all the flags have been set. The troubled-cells flagged for each variable is pooled together.

[back to table of contents](#markdown-header-table-of-contents)
