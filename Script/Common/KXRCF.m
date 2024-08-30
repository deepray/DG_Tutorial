function [ind] = KXRCF(u,vel,Mesh,KXRCF_M)


eps = 1e-10;

% uint = Mesh.hK.*sum(u(:,2:end-1).*(Mesh.M*u(:,2:end-1)))/2;
% uint = sqrt(uint);
uavg = Mesh.AVG1D*abs(u(:,2:end-1));
uL = abs(u(1,2:end-1)-u(end,1:end-2));
uR = abs(u(end,2:end-1)-u(1,3:end));
numerator = (uL.*(vel>0) + uR.*(vel<=0)).*(uavg>=eps);
denominator = (Mesh.hK.^((Mesh.N + 1)/2)).*uavg.*(uavg>=eps)+1.*(uavg<eps);


T = numerator./denominator;

% figure(100)
% %semilogy(max(T,1e-12))
% semilogy(max(uL,1e-12))
% %semilogy(max(numerator,1e-12),'-')
% %semilogy(max(denominator,1e-12),'--')


% T = numerator./denominator;
ind = find(T>KXRCF_M);

end