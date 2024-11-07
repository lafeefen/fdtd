function [kdx,kdy,kdz,bx,by,bz,axd,ayd,azd,bx1,by1,bz1,axd1,ayd1,azd1] = CPML_Region(Dimensions,Sep,CPML,dx,dy,dz,dt,Nx,Ny,Nz)
    fundamentalConstantsUnits;
    bx1=0;by1=0;bz1=0;
    ax1=0;ay1=0;az1=0;
    m=3;                % Scaling Order 2 to 6       | 0.85  for UPML
    ma=1;               % Scaling order ~1
    d = CPML*dy;
    RefCo = exp(-16);   % ReflectionCoefficient~exp(-16) | exp(-6) for UPML|
    SigMax = 0.75*((m+1)/(150*pi*dx));%-log(RefCo)*(m+1)/(2*Nu0*d);
    kxMax = 4;          % 7 to 20                      | 4.8  for UPML
    AlpMax = 0.000009;      % 0.15 to 0.3
    %==============| Fill the CPML regions (polynomial grading)|===============
    x=0*dy;
    kX = zeros(1,CPML);
    SigX = zeros(1,CPML);
    AlpX = zeros(1,CPML);
    for is = 1:CPML
        SigX(is) = SigMax*((x/d)^m);
        kX(is) = 1+(kxMax-1)*(x/d)^m;
        AlpX(is) = AlpMax*((d-x)/d)^ma;
        x=x+dy;
    end
    bX = exp(-((SigX./kX)+AlpX)*(dt/Ep0));
    aX = (SigX./(kX.*(SigX+(kX.*AlpX)))).*(bX-1);
    %{
    if(Dimensions==1)
        % {
        kx = 1; ky = 1; kz(1:Nz) = 1;
        bx = 1;by = 1;
        bz(1:Nz) = exp(-(AlpMax)*(dt/Ep0));
        az = (bz-1);
        for ii = 1:CPML
            bz(CPML+1-ii) = bX(ii);
            az(CPML+1-ii) = aX(ii);
            kz(CPML+1-ii) = kX(ii);
            bz(Nz-CPML+ii)= bX(ii);
            az(Nz-CPML+ii)= aX(ii);
            kz(Nz-CPML+ii)= kX(ii);
        end
        axd=0;ayd=0;
        azd = az/dy;
    end
    if(Dimensions==2)
        % {
        kx(1:Nx) = 1; ky(1:Ny) = 1; kz = 1;
        bx(1:Nx) = exp(-(AlpMax)*(dt/Ep0));
        by(1:Ny) = exp(-(AlpMax)*(dt/Ep0));
        bz = 1; ax = (bx-1); ay = (by-1);
        for ii = 1:CPML
            bx(CPML+1-ii) = bX(ii);
            by(CPML+1-ii) = bX(ii);
            ax(CPML+1-ii) = aX(ii);
            ay(CPML+1-ii) = aX(ii);
            kx(CPML+1-ii) = kX(ii);
            ky(CPML+1-ii) = kX(ii);
            bx(Nx-CPML+ii)= bX(ii);
            by(Ny-CPML+ii)= bX(ii);
            ax(Nx-CPML+ii)= aX(ii);
            ay(Ny-CPML+ii)= aX(ii);
            kx(Nz-CPML+ii)= kX(ii);
            ky(Nz-CPML+ii)= kX(ii);
        end
        axd = ax/dy;
        ayd = ay/dy;
        azd = 0;
    end
    if(Dimensions==3)
    %}
    % {
    kx(1:Nx,1) = 1; ky(1,1:Ny) = 1; kz(1,1,1:Nz) = 1;
    if Sep == 1
        bx(1:Nx+1,1) = 1;%exp(-(AlpMax)*(dt/Ep0));
        by(1,1:Ny+1) = 1;%exp(-(AlpMax)*(dt/Ep0));
        bz(1,1,1:Nz+1) = 1;%exp(-(AlpMax)*(dt/Ep0));
        ax = zeros(Nx+1,1); ay = zeros(1,Ny+1); az = zeros(1,1,Nz+1);
        for ii = 1:CPML
            bx(CPML+1-ii,1) = bX(ii);  by(1,CPML+1-ii) = bX(ii);  bz(1,1,CPML+1-ii) = bX(ii);
            ax(CPML+1-ii,1) = aX(ii);  ay(1,CPML+1-ii) = aX(ii);  az(1,1,CPML+1-ii) = aX(ii);
            kx(CPML+1-ii,1) = kX(ii);  ky(1,CPML+1-ii) = kX(ii);  kz(1,1,CPML+1-ii) = kX(ii);
            bx(Nx-CPML+ii,1) = bX(ii); by(1,Ny-CPML+ii) = bX(ii); bz(1,1,Nz-CPML+ii) = bX(ii);
            ax(Nx-CPML+ii,1) = aX(ii); ay(1,Ny-CPML+ii) = aX(ii); az(1,1,Nz-CPML+ii) = aX(ii);
            kx(Nx-CPML+ii,1) = kX(ii); ky(1,Ny-CPML+ii) = kX(ii); kz(1,1,Nz-CPML+ii) = kX(ii);
        end
    elseif Sep == 2
        %{
        bx1(1:Nx+1,1) = 1;%exp(-(AlpMax)*(dt/Ep0));
        by1(1,1:Ny+1) = 1;%exp(-(AlpMax)*(dt/Ep0));
        bz1(1,1,1:Nz+1) = 1;%exp(-(AlpMax)*(dt/Ep0));
        ax1 = zeros(Nx+1,1); ay1 = zeros(1,Ny+1); az1 = zeros(1,1,Nz+1);
        for ii = 1:CPML
            bx1(CPML+1-ii,1) = bX(ii);  by1(1,CPML+1-ii) = bX(ii);  bz1(1,1,CPML+1-ii) = bX(ii);
            ax1(CPML+1-ii,1) = aX(ii);  ay1(1,CPML+1-ii) = aX(ii);  az1(1,1,CPML+1-ii) = aX(ii);
            %kx1(CPML+1-ii,1) = kX(ii);  ky1(1,CPML+1-ii) = kX(ii);  kz1(1,1,CPML+1-ii) = kX(ii);
            bx1(Nx-CPML+ii,1) = bX(ii); by1(1,Ny-CPML+ii) = bX(ii); bz1(1,1,Nz-CPML+ii) = bX(ii);
            ax1(Nx-CPML+ii,1) = aX(ii); ay1(1,Ny-CPML+ii) = aX(ii); az1(1,1,Nz-CPML+ii) = aX(ii);
            %kx1(Nx-CPML+ii,1) = kX(ii); ky1(1,Ny-CPML+ii) = kX(ii); kz1(1,1,Nz-CPML+ii) = kX(ii);
        end
        %}
        bx(1:2*CPML,1) = 1;%exp(-(AlpMax)*(dt/Ep0));
        by(1,1:2*CPML) = 1;%exp(-(AlpMax)*(dt/Ep0));
        bz(1,1,1:2*CPML) = 1;%exp(-(AlpMax)*(dt/Ep0));
        ax = zeros(2*CPML,1); ay = zeros(1,2*CPML); az = zeros(1,1,2*CPML);
        for ii = 1:CPML
            bx(CPML+1-ii,1) = bX(ii);  by(1,CPML+1-ii) = bX(ii);  bz(1,1,CPML+1-ii) = bX(ii);
            ax(CPML+1-ii,1) = aX(ii);  ay(1,CPML+1-ii) = aX(ii);  az(1,1,CPML+1-ii) = aX(ii);
            kx(CPML+1-ii,1) = kX(ii);  ky(1,CPML+1-ii) = kX(ii);  kz(1,1,CPML+1-ii) = kX(ii);
            bx(CPML+ii,1) = bX(ii);    by(1,CPML+ii) = bX(ii);    bz(1,1,CPML+ii) = bX(ii);
            ax(CPML+ii,1) = aX(ii);    ay(1,CPML+ii) = aX(ii);    az(1,1,CPML+ii) = aX(ii);
            kx(Nx-CPML+ii,1) = kX(ii); ky(1,Ny-CPML+ii) = kX(ii); kz(1,1,Nz-CPML+ii) = kX(ii);
        end
    end
    kdx = dt./(kx.*dx);    kdy = dt./(ky.*dy);    kdz = dt./(kz.*dz);
    axd1 = ax1/dy;           ayd1 = ay1/dy;           azd1 = az1/dy;
    axd = ax/dy;           ayd = ay/dy;           azd = az/dy;
    %figure(4);plot(kx);figure(5);plot(kdx);
end