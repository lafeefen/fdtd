function [sx1,sy1,sz1,sx2,sy2,sz2,sxy1,sxy2,syz1,syz2,szx1,szx2] = UPML_Region_Sep(CPML,dx,dy,dz,dt,Nx,Ny,Nz,e0,u0)
[kx,ky,kz]=deal(1.0,1.0,1.0); [sx,sy,sz]=deal(0.0,0.0,0.0); tepdt = 2*e0*dt;
s1x(1:Nx+1,1) = (2*e0*kx - sx*dt);   s1y(1,1:Ny+1) = (2*e0*ky - sy*dt);   s1z(1,1,1:Nz+1) = (2*e0*kz - sz*dt);
s2x(1:Nx+1,1) = (2*e0*kx + sx*dt);   s2y(1,1:Ny+1) = (2*e0*ky + sy*dt);   s2z(1,1,1:Nz+1) = (2*e0*kz + sz*dt);
Nu0 = sqrt(u0/e0);                        % impedance of free space
RefCo = 1.0e-7;                           % reflectionCoefficient   ReflectionCoefficient~exp(-16)(CPML) | exp(-6) for UPML
m = 3;                                    % 2 <= m <= 6             Scaling Order 2 to 6(CPML)           | 0.85 (UPML)
d = CPML * dx;                            % width of PML region (in mm)
kxMax = 7;                                % 7 to 20                                                      | 4.8  (UPML)
SigMax = -log(RefCo)*(m+1.0)/(2.0*Nu0*d); % Maximum electrical Conductivity
SigBndryFactor = SigMax/(dx*(d^m)*(m+1));
%gradK = 1;
x = 0;
for ii = 1:CPML
    x1 = (x-0.5)*dx;
    x2 = (x+0.5)*dx;
    %gradSig = SigMax*(((x*dx)/d)^m);
    gradK = 1+(kxMax-1)*((x*dx)/d)^m;
    t1 = 2*e0*gradK;
    if ii == 1
        gradSig = SigBndryFactor*(x1^(m+1));
    else
        gradSig = SigBndryFactor*((x2^(m+1))-(x1^(m+1)));
    end
    t2 = gradSig*dt;
    s1x(CPML+1-ii) = t1 - t2;     s2x(CPML+1-ii) = t1 + t2;
    s1y(CPML+1-ii) = t1 - t2;     s2y(CPML+1-ii) = t1 + t2;
    s1z(CPML+1-ii) = t1 - t2;     s2z(CPML+1-ii) = t1 + t2;
    s1x(end-CPML+ii) = t1 - t2;   s2x(end-CPML+ii) = t1 + t2;
    s1y(end-CPML+ii) = t1 - t2;   s2y(end-CPML+ii) = t1 + t2;
    s1z(end-CPML+ii) = t1 - t2;   s2z(end-CPML+ii) = t1 + t2;
    x=x+1;
end

s2x_1 = 1./s2x; s2y_1 = 1./s2y; s2z_1 = 1./s2z;
[sx1,sy1,sz1,sx2,sy2,sz2] = deal( s1x./s2x, s1y./s2y, s1z./s2z, tepdt./(s2x*dx), tepdt./(s2y*dy), tepdt./(s2z*dz) );
[sxy1,sxy2] = deal( s2x*s2y_1,s1x*s2y_1 );
syz1(1,1:Ny+1,1:Nz+1) = s2y(1,Ny/2)*s2z_1(1,1,Nz/2); syz2(1,1:Ny+1,1:Nz+1) = s1y(1,Ny/2)*s2z_1(1,1,Nz/2);
szx1(1:Nx+1,1,1:Nz+1) = s2z(1,1,Nz/2)*s2x_1(Nx/2,1); szx2(1:Nx+1,1,1:Nz+1) = s1z(1,1,Nz/2)*s2x_1(Nx/2,1);
for ij=1:Ny+1
    for ik = 1:Nz+1
        syz1(1,ij,ik) = s2y(1,ij)*s2z_1(1,1,ik); syz2(1,ij,ik) = s1y(1,ij)*s2z_1(1,1,ik); % UPML Coeffcients
    end
end
for ii=1:Nx+1
    for ik = 1:Nz+1
        szx1(ii,1,ik) = s2z(1,1,ik)*s2x_1(ii,1); szx2(ii,1,ik) = s1z(1,1,ik)*s2x_1(ii,1); % UPML Coeffcients
    end
end
%{
%function [s1x,s1y,s1z,s2x,s2y,s2z] = UPML_Region_Sep(Dimensions,CPML,dz,dt,Nx,Ny,Nz)
    fundamentalConstantsUnits;
    % {
    m = 0.85;          % Scaling Order 2 to 6           | 0.85  for UPML
    d = CPML*dz;
    RefCo = exp(-6);   % ReflectionCoefficient~exp(-16) | exp(-6) for UPML
    SigMax = -log(RefCo)*(m+1)/(2*Nu0*d);
    kxMax = 4.8;       % 7 to 20                        | 4.8  for UPML
    %==============| Fill the CPML regions (polynomial grading)|===============
    x=0*dz;
    kX = zeros(1,CPML);
    SigX = zeros(1,CPML);
    for is = 1:CPML
        SigX(1,is) = SigMax*((x/d)^m);
        kX(1,is) = 1+(kxMax-1)*(x/d)^m;
        x=x+dz;
    end
    s1 = (kX/dt) - (SigX/tep);
    s2 = (kX/dt) + (SigX/tep);
    % }
    % {
    if(Dimensions==3)
        sxm(1,1:Nx) = 2*Ep0;          sym(1,1:Ny) = 2*Ep0;     szm(1,1:Nz) = 2*Ep0;
        sxp(1,1:Nx) = 2*Ep0;          syp(1,1:Ny) = 2*Ep0;     szp(1,1:Nz) = 2*Ep0;
        sxm(1,1:CPML) = flip(s1,1);   sxp(1,1:CPML) = flip(s2,1);
        sxm(1,(Nx-CPML+1):Nx) = s1;   sxp(1,(Nx-CPML+1):Nx) = s2;  sxm(1,Nx+1) = sxm(1,Nx);  sxp(1,Nx+1) = sxp(1,Nx);
        sym(1,1:CPML) = flip(s1,1);   syp(1,1:CPML) = flip(s2,1);
        sym(1,(Ny-CPML+1):Ny) = s1;   syp(1,(Ny-CPML+1):Ny) = s2;  sym(1,Ny+1) = sym(1,Ny);  syp(1,Ny+1) = syp(1,Ny);
        szm(1,1:CPML) = flip(s1,1);   szp(1,1:CPML) = flip(s2,1);
        szm(1,(Nz-CPML+1):Nz) = s1;   szp(1,(Nz-CPML+1):Nz) = s2;  szm(1,Nz+1) = szm(1,Nz);  szp(1,Nz+1) = szp(1,Nz);
        sxp_1 = 1./sxp;               syp_1 = 1./syp;              szp_1 = 1./szp;

        sx1 = sxm.*sxp_1;             sx2 = 2*Ep0*dt.*sxp_1/dz;
        sy1 = sym.*syp_1;             sy2 = 2*Ep0*dt.*syp_1/dz;
        sz1 = szm.*szp_1;             sz2 = 2*Ep0*dt.*szp_1/dz;
        sxy1 = (sxp.')*syp_1;         sxy2 = (sxm.')*syp_1;
        syz1 = (syp.')*szp_1;         syz2 = (sym.')*szp_1;
        szx1 = (szp.')*sxp_1;         szx2 = (szm.')*sxp_1;
    end
    %}
end