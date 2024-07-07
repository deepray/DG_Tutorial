function [u] = Scalar1D_Visc(u,Problem,Mesh,Viscosity,Output)

% Purpose  : Integrate 1D Scalar equation until
%            FinalTime starting with
%            initial condition u in the domain [xL,xR].

% Set exact flux and Jacobian
[flux,dflux] = Set_scalar_flux1D(Problem.model);

time = 0;
iter = 0;
dxmin = min(Mesh.hK);
dt = 0;

fid = fopen(strcat(Output.fname_base,'_visc.dat'),'w');


figure(1)
plot(Mesh.x(:),u(:),'b-','LineWidth',2)
xlabel('x')
ylabel('u')
xlim([Mesh.bnd_l,Mesh.bnd_r])
title(['time = ',num2str(time)])

% Initialize solution at previous time step
u_tmp=zeros(length(u(:)),3);

% outer time step loop
while(time<Problem.FinalTime)

    %Compute artificial viscosity
    uold = reshape(u_tmp(:,1),Mesh.Np,Mesh.K);
    mu_piece = Scalar1D_viscosity(u, uold, dflux, Viscosity, Problem, Mesh, dt, iter);
    mu_piece = max(mu_piece,0);
    mu_vals  = Scalar1D_smooth_viscosity(mu_piece,Mesh.x); 
    mu_vals  = max(mu_vals,0);
    maxvisc  = max(abs(mu_vals(:)));
    
    if(mod(iter,Output.ind_iter) == 0 || time+dt >= Problem.FinalTime)
        Visc_write1D(fid,time,mu_vals);
    end
    
    speed = max(max(abs(dflux(u))));
    dt = Problem.CFL*1/(speed*Mesh.N^2/dxmin+maxvisc*Mesh.N^4/dxmin^2);
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % 3rd order SSP Runge-Kutta
    
    % SSP RK Stage 1.
    rhsu  = ScalarRHS1D_viscweak(u,flux,dflux,mu_vals,Problem.bc_cond,Mesh);
    u1  = u  + dt*rhsu;
    
    % SSP RK Stage 2.
    rhsu  = ScalarRHS1D_viscweak(u1,flux,dflux,mu_vals,Problem.bc_cond,Mesh);
    u2   = (3*u  + u1  + dt*rhsu )/4;
    
 
    % SSP RK Stage 3.
    rhsu  = ScalarRHS1D_viscweak(u2,flux,dflux,mu_vals,Problem.bc_cond,Mesh);
    u  = (u  + 2*u2  + 2*dt*rhsu )/3;
    
    % Increment time and adapt timestep
    time = time+dt;

    % Increment saved variables (needed for EV)
    u_tmp(:,3)=u(:); u_tmp(:,1)=u_tmp(:,2); u_tmp(:,2)=u_tmp(:,3);
    
    if(mod(iter,Output.plot_iter) == 0 || time >= Problem.FinalTime)
        
        figure(1)
        plot(Mesh.x(:),u(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('u')
        title(['time = ',num2str(time)])
        pause(0.1)
    end
    
    iter = iter + 1;
    
end

fclose(fid);

return
