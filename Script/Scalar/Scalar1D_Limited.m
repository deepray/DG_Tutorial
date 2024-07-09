function [u] = Scalar1D_Limited(u,Problem,Mesh,Limit,Net,Output)

% Purpose  : Integrate 1D Scalar equation until
%            FinalTime starting with
%            initial condition u in the domain [xL,xR].

% Set exact flux and Jacobian
[flux,dflux] = Set_scalar_flux1D(Problem.model);

time = 0;
iter = 0;
xcen = mean(Mesh.x,1);
dxmin = min(Mesh.hK);

fid = fopen(strcat(Output.fname_base,'_tcells.dat'),'w');


% Limit initial solution
ind0  = Scalar1D_Tcells(u,Problem.bc_cond,Mesh,Limit,Net);
u     = Scalar1D_limit(u,ind0,Problem.bc_cond,Limit.Limiter,Mesh);


Tcell_write1D(fid,time,xcen(ind0));
figure(1)
subplot(2,1,1)
plot(Mesh.x(:),u(:),'b-','LineWidth',2)
xlabel('x')
ylabel('u')
xlim([Mesh.bnd_l,Mesh.bnd_r])
title(['time = ',num2str(time)])

subplot(2,1,2)
plot(xcen(ind0),ones(1,length(ind0))*time,'r.')
xlabel('x')
ylabel('t')
xlim([Mesh.bnd_l,Mesh.bnd_r])
ylim([0,Problem.FinalTime])
hold all


% outer time step loop
while(time<Problem.FinalTime)
    
    speed = max(max(abs(dflux(u))));
    dt = Problem.CFL* min(dxmin/abs(speed))/(Mesh.N^2);
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % 3rd order SSP Runge-Kutta
    
    % SSP RK Stage 1.
    rhsu  = ScalarRHS1D_weak(u,flux,dflux,Problem.bc_cond,Mesh);
    u1  = u  + dt*rhsu;
    
    % Limit fields
    ind1  = Scalar1D_Tcells(u1,Problem.bc_cond,Mesh,Limit,Net);
    u1    = Scalar1D_limit(u1,ind1,Problem.bc_cond,Limit.Limiter,Mesh);
    
    
    % SSP RK Stage 2.
    rhsu  = ScalarRHS1D_weak(u1,flux,dflux,Problem.bc_cond,Mesh);
    u2   = (3*u  + u1  + dt*rhsu )/4;
    
    % Limit fields
    ind2  = Scalar1D_Tcells(u2,Problem.bc_cond,Mesh,Limit,Net);
    u2    = Scalar1D_limit(u2,ind2,Problem.bc_cond,Limit.Limiter,Mesh);
    
    % SSP RK Stage 3.
    rhsu  = ScalarRHS1D_weak(u2,flux,dflux,Problem.bc_cond,Mesh);
    u  = (u  + 2*u2  + 2*dt*rhsu )/3;
    
    % Limit solution
    ind3  = Scalar1D_Tcells(u,Problem.bc_cond,Mesh,Limit,Net);
    u     = Scalar1D_limit(u,ind3,Problem.bc_cond,Limit.Limiter,Mesh);
    
    inda = unique([ind1,ind2,ind3]);
    Tcell_write1D(fid,time+dt,xcen(inda));

    
    % Increment time and adapt timestep
    time = time+dt;
    
    if(mod(iter,Output.plot_iter) == 0 || time >= Problem.FinalTime)
        
        figure(1)
        subplot(2,1,1)
        plot(Mesh.x(:),u(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('u')
        
        subplot(2,1,2)
        plot(xcen(inda),ones(1,length(inda))*time,'r.') 
        pause(0.1)
    end
    
    iter = iter + 1;
    
end

fclose(fid);

return
