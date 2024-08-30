% Check parameters
ScalarCheckParam1D;

% Display paramaters
ScalarStartDisp1D;

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
u = Problem.u_IC(Mesh.x);

% Creating file name base for saving solution
Output.fname_base = Scalar_fnamebase1D(Problem,Mesh.N,Mesh.K,Limit,Viscosity,Mesh.mesh_pert);

% Solve Problem
tic
fprintf('... starting main solve\n')
if(~strcmp(Limiter,'NONE') |  strcmp(Visc_model,'NONE'))
    [u] = Scalar1D_Limited(u,Problem,Mesh,Limit,Net,Output);
else
    [u] = Scalar1D_Visc(u,Problem,Mesh,Viscosity,Output);
end
toc

% Save final solution
Save_scalar_soln1D(u,Mesh.x,Output.fname_base);

fprintf('... generating and saving plots in directory OUTPUT\n')
if(Output.ref_avail)
    PlotScalar1D(Output.fname_base,Mesh.x,[Mesh.bnd_l,Mesh.bnd_r],Output.var_ran,[0,Problem.FinalTime],Output.ind_iter,true,Output.ref_fname); 
else
    PlotScalar1D(Output.fname_base,Mesh.x,[Mesh.bnd_l,Mesh.bnd_r],Output.var_ran,[0,Problem.FinalTime],Output.ind_iter,false);
end


% Clean up
fprintf('... cleaning up\n')
CleanUp1D;

fprintf('------------ Solver has finished -------------\n')


