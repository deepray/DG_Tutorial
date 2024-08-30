function ind = Scalar1D_Tcells(u,bc_cond,Mesh,Limit,Net,dflux)

u_ext = Apply_BC1D(u,bc_cond);

vel   = Mesh.AVG1D*dflux(u);

ind   = Find_Tcells1D(u_ext,vel,Mesh,Limit,Net);

return

