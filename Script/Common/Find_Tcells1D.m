function ind = Find_Tcells1D(u,Mesh,Limit,Net)

% Purpose: find all the troubled-cells for variable u
% NOTE: u must include ghost cell values

% Compute cell averages
v = Mesh.AVG1D*u;

eps0=1.0e-8;

% find end values of each element (excluding ghost cells)
ue1 = u(1,2:end-1); ue2 = u(end,2:end-1);

% find cell averages 
vk = v(2:Mesh.K+1); vkm1 = v(1:Mesh.K); vkp1 = v(3:Mesh.K+2);

% Find elements in need of limiting
if(strcmp(Limit.Indicator,'NONE'))
    ind = [];
elseif(strcmp(Limit.Indicator,'ALL'))
    ind = 1:Mesh.K;    
elseif(strcmp(Limit.Indicator,'MINMOD'))
    ve1 = vk - minmod([(vk-ue1);vk-vkm1;vkp1-vk]);
    ve2 = vk + minmod([(ue2-vk);vk-vkm1;vkp1-vk]);
    ind = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);
elseif(strcmp(Limit.Indicator,'TVB'))
    ve1 = vk - minmodB([(vk-ue1);vk-vkm1;vkp1-vk],Limit.TVBM,Mesh.x(end,:)-Mesh.x(1,:));
    ve2 = vk + minmodB([(ue2-vk);vk-vkm1;vkp1-vk],Limit.TVBM,Mesh.x(end,:)-Mesh.x(1,:));
    ind = find(abs(ve1-ue1)>eps0 | abs(ve2-ue2)>eps0);
elseif(strcmp(Limit.Indicator,'NN'))
    tcell = ind_MLP1D([vkm1;vk;vkp1;ue1;ue2],Net);
    ind = find(tcell==1); 
else
    error('Indicator %s not available!!',Limit.Indicator)
end


return
