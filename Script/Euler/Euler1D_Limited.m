function q = Euler1D_Limited(q,Problem,Mesh,Limit,Net,Output)

% Purpose  : Integrate 1D Euler equations until FinalTime

time = 0;

dxmin = min(abs(Mesh.x(1,:)-Mesh.x(2,:)));
iter = 0;
xcen = mean(Mesh.x,1);

% Limit initial solution
fid = fopen(strcat(Output.fname_base,'_tcells.dat'),'w');

ind0 = Tcells_Euler_type1D(q,Problem,Mesh,Limit,Net);
q    = SlopeLimit_Euler_type1D(q,ind0,Problem,Limit,Mesh);

Tcell_write1D(fid,time,xcen(ind0));

density = q(:,:,1);
vel     = q(:,:,2)./q(:,:,1);
pre     = (Problem.gas_gamma-1)*(q(:,:,3) - 0.5*q(:,:,2).^2./q(:,:,1));

figure(1)
subplot(2,3,1)
plot(Mesh.x(:),density(:),'b-','LineWidth',2)
xlabel('x')
ylabel('Density')
title(['time = ',num2str(time)])

subplot(2,3,2)
plot(Mesh.x(:),vel(:),'b-','LineWidth',2)
xlabel('x')
ylabel('Velocity')

subplot(2,3,3)
plot(Mesh.x(:),pre(:),'b-','LineWidth',2)
xlabel('x')
ylabel('Pressure')

subplot(2,3,4)
plot(xcen(ind0),ones(1,length(ind0))*time,'r.')
xlabel('x')
ylabel('t')
xlim([Mesh.bnd_l Mesh.bnd_r])
ylim([0 Problem.FinalTime])
hold all

subplot(2,3,5)
plot(xcen(ind0),ones(1,length(ind0))*time,'r.')
xlabel('x')
ylabel('t')
xlim([Mesh.bnd_l Mesh.bnd_r])
ylim([0 Problem.FinalTime])
hold all

subplot(2,3,6)
plot(xcen(ind0),ones(1,length(ind0))*time,'r.')
xlabel('x')
ylabel('t')
xlim([Mesh.bnd_l Mesh.bnd_r])
ylim([0 Problem.FinalTime])
hold all

% outer time step loop
while(time<Problem.FinalTime)
    
    lambda = sqrt(Problem.gas_gamma*pre./q(:,:,1)) + abs(q(:,:,2)./q(:,:,1));
    dt     = Problem.CFL*min(min(dxmin./(lambda)));
    
    if(time+dt>Problem.FinalTime)
        dt = Problem.FinalTime-time;
    end
    
    % 3rd order SSP Runge-Kutta
    
    % SSP RK Stage 1.
    [rhsq]  = EulerRHS1D_weak(q, Problem.gas_gamma, Problem.gas_const,Problem.bc_cond,Mesh);
    q1      = q + dt*rhsq;
    
    ind1 = Tcells_Euler_type1D(q1,Problem,Mesh,Limit,Net);
    q1   = SlopeLimit_Euler_type1D(q1,ind1,Problem,Limit,Mesh);
    
    
    pre = (Problem.gas_gamma-1)*(q1(:,:,3) - 0.5*q1(:,:,2).^2./q1(:,:,1));
    if( min(min(real(q1(:,:,1)))) <= 0.0 || min(min(real(pre))) <= 0.0)
        error('Positivity loss!!');
    end
    
    
    % SSP RK Stage 2.
    [rhsq]  = EulerRHS1D_weak(q1, Problem.gas_gamma, Problem.gas_const,Problem.bc_cond,Mesh);
    q2      = (3*q + (q1 + dt*rhsq))/4.0;
    
    ind2 = Tcells_Euler_type1D(q2,Problem,Mesh,Limit,Net);
    q2   = SlopeLimit_Euler_type1D(q2,ind2,Problem,Limit,Mesh);
    
    
    pre    = (Problem.gas_gamma-1)*(q2(:,:,3) - 0.5*q2(:,:,2).^2./q2(:,:,1));
    if( min(min(real(q2(:,:,1)))) <= 0.0 || min(min(real(pre))) <= 0.0)
        error('Positivity loss!!');
    end
    
    
    % SSP RK Stage 3.
    [rhsq]  = EulerRHS1D_weak(q2,Problem.gas_gamma, Problem.gas_const,Problem.bc_cond,Mesh);
    q       = (q + 2*(q2 + dt*rhsq))/3.0;
    
    ind3 = Tcells_Euler_type1D(q,Problem,Mesh,Limit,Net);
    q    = SlopeLimit_Euler_type1D(q,ind3,Problem,Limit,Mesh);
    

    inda = unique([ind1,ind2,ind3]);
    Tcell_write1D(fid,time+dt,xcen(inda));
    
    
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
        subplot(2,3,1)
        plot(Mesh.x(:),density(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Density')
        title(['time = ',num2str(time)])
        
        subplot(2,3,2)
        plot(Mesh.x(:),vel(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Velocity')
        
        subplot(2,3,3)
        plot(Mesh.x(:),pre(:),'b-','LineWidth',2)
        xlabel('x')
        ylabel('Pressure')
        
        subplot(2,3,4)
        plot(xcen(inda),ones(1,length(inda))*time,'r.')
        xlabel('x')
        ylabel('t')
        xlim([Mesh.bnd_l Mesh.bnd_r])
        ylim([0 Problem.FinalTime])
        hold all
        
        subplot(2,3,5)
        plot(xcen(inda),ones(1,length(inda))*time,'r.')
        xlabel('x')
        ylabel('t')
        xlim([Mesh.bnd_l Mesh.bnd_r])
        ylim([0 Problem.FinalTime])
        hold all
        
        subplot(2,3,6)
        plot(xcen(inda),ones(1,length(inda))*time,'r.')
        xlabel('x')
        ylabel('t')
        xlim([Mesh.bnd_l Mesh.bnd_r])
        ylim([0 Problem.FinalTime])
        hold all
        
        pause(.01)
        
    end  
    iter = iter + 1;
    
end

fclose(fid);

return
