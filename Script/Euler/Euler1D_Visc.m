function q = Euler1D_Visc(q,Problem,Mesh,Viscosity,Output)

% Purpose  : Integrate 1D Euler equations until FinalTime

time = 0;
dt   = 0;

dxmin = min(Mesh.hK);
iter = 0;

% Limit initial solution
fid = fopen(strcat(Output.fname_base,'_visc.dat'),'w');

density = q(:,:,1);
vel     = q(:,:,2)./q(:,:,1);
pre     = (Problem.gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1));

figure(1)
subplot(1,3,1)
plot(Mesh.x(:),density(:),'b-','LineWidth',2)
xlabel('x')
ylabel('Density')
title(['time = ',num2str(time)])

subplot(1,3,2)
plot(Mesh.x(:),vel(:),'b-','LineWidth',2)
xlabel('x')
ylabel('Velocity')

subplot(1,3,3)
plot(Mesh.x(:),pre(:),'b-','LineWidth',2)
xlabel('x')
ylabel('Pressure')

% Initialize solution at previous time step
q_tmp=zeros(length(q(:)),3);

% outer time step loop
while(time<Problem.FinalTime)
    
    %Compute artificial viscosity
    qold = reshape(q_tmp(:,1),Mesh.Np,Mesh.K,3);
    mu_piece = Euler1D_viscosity(q, qold, Viscosity, Problem, Mesh, dt, iter);
    mu_piece=max(mu_piece,0);
    mu_vals=Scalar1D_smooth_viscosity(mu_piece,Mesh.x); 
    mu_vals=max(mu_vals,0);
    maxvisc=max(abs(mu_vals(:)));

    if(mod(iter,Output.ind_iter) == 0 || time+dt >= Problem.FinalTime)
        Visc_write1D(fid,time,mu_vals);
    end

    %Apply same AV for all equations
    mu_vals=repmat(mu_vals,1,1,3);
    
    
    lambda = sqrt(Problem.gas_gamma*pre./q(:,:,1)) + abs(q(:,:,2)./q(:,:,1));
    dt = Problem.CFL*1/(max(lambda(:))*Mesh.N^2/dxmin+maxvisc*Mesh.N^4/dxmin^2);
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % 3rd order SSP Runge-Kutta
    
    % SSP RK Stage 1.
    [rhsq]  = EulerRHS1D_viscweak(q, Problem.gas_gamma, Problem.gas_const,mu_vals,Problem.bc_cond,Mesh);
    q1      = q + dt*rhsq;
        
    pre = (Problem.gas_gamma-1)*(q1(:,:,3) - 0.5*q1(:,:,2).^2./q1(:,:,1));
    if( min(min(real(q1(:,:,1)))) <= 0.0 || min(min(real(pre))) <= 0.0)
        error('Positivity loss!!');
    end
    
    
    % SSP RK Stage 2.
    [rhsq]  = EulerRHS1D_viscweak(q1, Problem.gas_gamma, Problem.gas_const,mu_vals,Problem.bc_cond,Mesh);
    q2      = (3*q + (q1 + dt*rhsq))/4.0;
        
    pre    = (Problem.gas_gamma-1)*(q2(:,:,3) - 0.5*q2(:,:,2).^2./q2(:,:,1));
    if( min(min(real(q2(:,:,1)))) <= 0.0 || min(min(real(pre))) <= 0.0)
        error('Positivity loss!!');
    end
    
    
    % SSP RK Stage 3.
    [rhsq]  = EulerRHS1D_viscweak(q2,Problem.gas_gamma, Problem.gas_const,mu_vals,Problem.bc_cond,Mesh);
    q       = (q + 2*(q2 + dt*rhsq))/3.0;
       
    pre    = (Problem.gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1));
    if( min(min(real(q(:,:,1)))) <= 0.0 || min(min(real(pre))) <= 0.0)
        error('Positivity loss!!');
    end
    
    % Increment time and adapt timestep
    time = time+dt;    
    
    if(mod(iter,Output.plot_iter) == 0 || time >= Problem.FinalTime)
        density = q(:,:,1);
        vel   = q(:,:,2)./q(:,:,1);
        figure(1)
        subplot(1,3,1)
        plot(Mesh.x(:),density(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Density')
        title(['time = ',num2str(time)])
        
        subplot(1,3,2)
        plot(Mesh.x(:),vel(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Velocity')
        
        subplot(1,3,3)
        plot(Mesh.x(:),pre(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Pressure')
        
        pause(.01)
        
    end  
    iter = iter + 1;
    
end

fclose(fid);

return
