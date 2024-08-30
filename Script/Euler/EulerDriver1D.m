% Check parameters
EulerCheckParam1D;

% Display paramaters
EulerStartDisp1D;

% Find relative path
Find_relative_path;

% Generate simple mesh
[Mesh.Nv, Mesh.VX, Mesh.hK] = MeshGen1D(Mesh.bnd_l,Mesh.bnd_r,Mesh.K, Mesh.mesh_pert);

% generate various matrix operators and maps
StartUp1D;

% Extract MLP weights, biases and other parameters
if(strcmp(Limit.Indicator,'NN'))
    Net = read_mlp_param1D(REL_PATH,Limit.NN_model);
else
    Net.avail = false;
end


% Generate mass matrix and initialize solution
rho        = Problem.rho_IC(Mesh.x);
vel        = Problem.vel_IC(Mesh.x);
pre        = Problem.pre_IC(Mesh.x);
mmt        = rho.*vel;
energy     = 0.5*rho.*vel.^2 + pre/(Problem.gas_gamma - 1);

% Creating vector of conserved variables
q = zeros(Mesh.Np,Mesh.K,3); q(:,:,1) = rho; q(:,:,2) = mmt;  q(:,:,3) = energy;

% Creating file name base for saving solution
Output.fname_base = Euler_fnamebase1D(Problem,Mesh.N,Mesh.K,Limit,Viscosity,Mesh.mesh_pert);


% Solve Problem
tic
fprintf('... starting main solve\n')
if(~strcmp(Limiter,'NONE'))
    q = Euler1D_Limited(q,Problem,Mesh,Limit,Net,Output);
else
    q = Euler1D_Visc(q,Problem,Mesh,Viscosity,Output);
end
toc


% Save final solution
Save_Euler_soln1D(q,Mesh.x,Problem.gas_gamma,Output.fname_base);

%%
fprintf('... generating and saving plots in directory OUTPUT\n')
if(Output.ref_avail)
    PlotEuler1D(Output.fname_base,Mesh.x,[Mesh.bnd_l,Mesh.bnd_r],Output.var_ran,[0,Problem.FinalTime],Output.ind_iter,true,Output.ref_fname); 
else
    PlotEulerE1D(Output.fname_base,Mesh.x,[Mesh.bnd_l,Mesh.bnd_r],Output.var_ran,[0,Problem.FinalTime],Output.ind_iter,false);
end


% Clean up
fprintf('... cleaning up\n')
CleanUp1D;

fprintf('------------ Solver has finished -------------\n')

