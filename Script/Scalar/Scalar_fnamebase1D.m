function fname = Scalar_fnamebase1D(Problem,N,K,Limit,Viscosity,mesh_pert)

% Globals1D_DG;
% Globals1D_MLP;

mkdir('OUTPUT');

fname = sprintf('OUTPUT/%s1D_%s_P%d_N%d',Problem.model,Problem.test_name,N,K);
if(strcmp(Limit.Indicator,'NONE'))
    fname = sprintf('%s_IND_%s',fname,Limit.Indicator);
elseif(strcmp(Limit.Indicator,'ALL'))
    fname = sprintf('%s_IND_%s',fname,Limit.Indicator);    
elseif(strcmp(Limit.Indicator,'MINMOD'))
    fname = sprintf('%s_IND_%s',fname,Limit.Indicator);
elseif(strcmp(Limit.Indicator,'TVB'))
    fname = sprintf('%s_IND_%s_%d',fname,Limit.Indicator,Limit.TVBM);
elseif(strcmp(Limit.Indicator,'NN'))
    fname = sprintf('%s_IND_%s',fname,Limit.Indicator);
elseif(strcmp(Limit.Indicator,'FuShu'))
    fname = sprintf('%s_IND_FS',fname);    
else
    error('Indicator %s not available',Limit.Indicator);
end
fname = sprintf('%s_LIM_%s',fname,Limit.Limiter);

if(strcmp(Viscosity.model,'NONE'))
    fname = sprintf('%s_VISC_%s',fname,Viscosity.model);
elseif(strcmp(Viscosity.model,'MDH'))
    fname = sprintf('%s_VISC_%s_%.3f_%.3f_%.3f',fname,Viscosity.model,Viscosity.c_A,Viscosity.c_k,Viscosity.c_max);
    %fname = sprintf('%s_VISC_%s',fname,Viscosity.model);
elseif(strcmp(Viscosity.model,'MDA'))
    fname = sprintf('%s_VISC_%s_%.3f',fname,Viscosity.model,Viscosity.c_max);
    %fname = sprintf('%s_VISC_%s',fname,Viscosity.model);
elseif(strcmp(Viscosity.model,'EV'))
    fname = sprintf('%s_VISC_%s_%.3f_%.3f',fname,Viscosity.model,Viscosity.c_E,Viscosity.c_max);
    %fname = sprintf('%s_VISC_%s',fname,Viscosity.model);
else
    error('Viscosity model %s not available',Viscosity.model);
end

if(mesh_pert ~= 0.0)
    fname = sprintf('%s_pert',fname);
end


return
