function PlotScalar1D(fname_base,x,x_ran,u_ran,t_ran,ind_iter,ref_avail,ref_fname)


soln = load(strcat(fname_base,'.dat'));

if(ref_avail) 
    ref_soln = load(ref_fname);
end

figure(10)
clf
if(ref_avail)
    plot(ref_soln(:,1),ref_soln(:,2),'k-','LineWidth',2) 
    hold all
end
plot(soln(:,1),soln(:,2),'r-','LineWidth',2)
if(ref_avail)
    legend('Reference','Numerical')
end
xlim(x_ran);
ylim(u_ran);
xlabel('x')
ylabel('u')
set(gca,'FontSize',20)
fname = sprintf('%s_soln.png',fname_base);
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
    while ischar(tline)
        data  = str2num(char(regexp(tline,', ','split')));
        t     = data(1);
        xc    = data(2:end);
        
        plot(xc,t*ones(1,length(xc)),'r.')   
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



