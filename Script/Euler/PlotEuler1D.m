function PlotEuler1D(fname_base,x,x_ran,var_ran,t_ran,ind_iter,ref_avail,ref_fname)


soln = load(strcat(fname_base,'.dat'));

if(ref_avail) 
    ref_soln = load(ref_fname);
end

figure(10)
clf

figure(11)
clf

figure(12)
clf
if(ref_avail)
    
    figure(10)
    plot(ref_soln(:,1),ref_soln(:,2),'k-','LineWidth',2) 
    hold all
    
    figure(11)
    plot(ref_soln(:,1),ref_soln(:,3),'k-','LineWidth',2) 
    hold all
    
    figure(12)
    plot(ref_soln(:,1),ref_soln(:,4),'k-','LineWidth',2) 
    hold all
end
figure(10)
plot(soln(:,1),soln(:,2),'r-','LineWidth',2)
figure(11)
plot(soln(:,1),soln(:,3),'r-','LineWidth',2)
figure(12)
plot(soln(:,1),soln(:,4),'r-','LineWidth',2)
if(ref_avail)
    figure(10)
    legend('Reference','Numerical')
    
    figure(11)
    legend('Reference','Numerical')
    
    figure(12)
    legend('Reference','Numerical')
end

figure(10)
xlim(x_ran);
ylim(var_ran(1,:));
xlabel('x')
ylabel('density')
set(gca,'FontSize',20)
fname = sprintf('%s_soln_density.png',fname_base);
print(fname,'-dpng')

figure(11)
xlim(x_ran);
ylim(var_ran(2,:));
xlabel('x')
ylabel('velocity')
set(gca,'FontSize',20)
fname = sprintf('%s_soln_vel.png',fname_base);
print(fname,'-dpng')

figure(12)
xlim(x_ran);
ylim(var_ran(3,:));
xlabel('x')
ylabel('pressure')
set(gca,'FontSize',20)
fname = sprintf('%s_soln_pre.png',fname_base);
print(fname,'-dpng')


% Plotting troubled cells if file exists
tcell_fname = strcat(fname_base,'_tcells.dat');
if(isfile(tcell_fname))
    fid = fopen(strcat(fname_base,'_tcells.dat'));
    
    figure(11)
    clf
    xlim(x_ran)
    ylim(t_ran)
    xlabel('x')
    ylabel('t')
    set(gca,'FontSize',20)
    hold all
    
    tline = fgetl(fid);
    iter = 0; 
    while ischar(tline)
        data  = str2num(char(regexp(tline,', ','split')));
        t     = data(1);
        xc    = data(2:end);
        
        if(mod(iter,ind_iter) == 0)
            plot(xc,t*ones(1,length(xc)),'r.')   
        end
        iter = iter + 1;
        tline = fgetl(fid);
    end
    fname = sprintf('%s_tcells.png',fname_base);
    print(fname,'-dpng')
    
    fclose(fid);
end

% Plotting viscosity if file exists
visc_fname = strcat(fname_base,'_visc.dat');
if(isfile(visc_fname))
    fid = fopen(strcat(fname_base,'_visc.dat'));
    tline = fgetl(fid);
    visc_= [];
    t = [];
    while ischar(tline)
        data  = str2num(char(regexp(tline,', ','split')));
        t     = [t; data(1)];
        visc_ = [visc_, data(2:end)];
        tline = fgetl(fid);
    end
    figure(11)
    clf
    visc_=log(max(visc_,1e-5))/log(10);
    plot_visc=pcolor(x(:),t,visc_'); shading interp; set(plot_visc, 'EdgeColor', 'none'); colorbar; colormap jet;
    xlabel({'$x$'},'interpreter','latex','FontSize',20); 
    ylabel({'$t$'},'interpreter','latex','FontSize',20);
    
    fname = sprintf('%s_visc.png',fname_base);
    print(fname,'-dpng')

    fclose(fid);
end



