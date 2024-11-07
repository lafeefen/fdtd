%% Start loop
tic
for t = 1:nt
    %% Magnetic Field Update
    Esrc_x = Ex_srcArray(t)*G_Distrtn;
    Esrc_y = Ey_srcArray(t)*G_Distrtn;
    ik = zSrc;
    DT(:,1:Ny,1:Nz) = Bx;
    CHyz(:,1:CPML,:) = by(1,1:CPML).*CHyz(:,1:CPML,:) + ayd(1,1:CPML).*diff(Ez(:,1:CPML+1,:),1,2);
    CHyz(:,CPML+1:end,:) = by(1,CPML+1:end).*CHyz(:,CPML+1:end,:) + ayd(1,CPML+1:end).*diff(Ez(:,Ny-CPML+1:Ny+1,:),1,2);
    CHzy(:,:,1:CPML) = bz(1,1,1:CPML).*CHzy(:,:,1:CPML) + azd(1,1,1:CPML).*diff(Ey(:,:,1:CPML+1),1,3);
    CHzy(:,:,CPML+1:end) = bz(1,1,CPML+1:end).*CHzy(:,:,CPML+1:end) + azd(1,1,CPML+1:end).*diff(Ey(:,:,Nz-CPML+1:Nz+1),1,3);
    Bx = Bx - kdy(1,1:Ny).*diff(Ez,1,2) + kdz(1,1,1:Nz).*diff(Ey,1,3);
    %Bx(:,:,ik) = Bx(:,:,ik) - kdz(1,1,ik).*Esrc_y(:,1:Ny);    % TFSF
    Bx(:,1:CPML,:) = Bx(:,1:CPML,:) - dt*CHyz(:,1:CPML,:);    Bx(:,Ny-CPML+1:Ny,:) = Bx(:,Ny-CPML+1:Ny,:) - dt*CHyz(:,CPML+1:end,:);
    Bx(:,:,1:CPML) = Bx(:,:,1:CPML) + dt*CHzy(:,:,1:CPML);    Bx(:,:,Nz-CPML+1:Nz) = Bx(:,:,Nz-CPML+1:Nz) + dt*CHzy(:,:,CPML+1:end);
    Hx = Hx + u_1.*( Bx - DT(:,1:Ny,1:Nz) );
    DT(1:Nx,:,1:Nz)= By;
    CHzx(:,:,1:CPML) = bz(1,1,1:CPML).*CHzx(:,:,1:CPML) + azd(1,1,1:CPML).*diff(Ex(:,:,1:CPML+1),1,3);
    CHzx(:,:,CPML+1:end) = bz(1,1,CPML+1:end).*CHzx(:,:,CPML+1:end) + azd(1,1,CPML+1:end).*diff(Ex(:,:,Nz-CPML+1:Nz+1),1,3);
    CHxz(1:CPML,:,:) = bx(1:CPML,1).*CHxz(1:CPML,:,:) + axd(1:CPML,1).*diff(Ez(1:CPML+1,:,:),1,1);
    CHxz(CPML+1:end,:,:) = bx(CPML+1:end,1).*CHxz(CPML+1:end,:,:) + axd(CPML+1:end,1).*diff(Ez(Nx-CPML+1:Nx+1,:,:),1,1);
    By = By - kdz(1,1,1:Nz).*diff(Ex,1,3) + kdx(1:Nx,1).*diff(Ez,1,1);
    %By(:,:,ik) = By(:,:,ik) + kdz(1,1,ik).*Esrc_y(1:Nx,:);   % TFSF
    By(:,:,1:CPML) = By(:,:,1:CPML) - dt*CHzx(:,:,1:CPML);   By(:,:,Nz-CPML+1:Nz) = By(:,:,Nz-CPML+1:Nz) - dt*CHzx(:,:,CPML+1:end);
    By(1:CPML,:,:) = By(1:CPML,:,:) + dt*CHxz(1:CPML,:,:);   By(Nx-CPML+1:Nx,:,:) = By(Nx-CPML+1:Nx,:,:) + dt*CHxz(CPML+1:end,:,:);
    Hy = Hy + u_1.*( By - DT(1:Nx,:,1:Nz) );
    DT(1:Nx,1:Ny,:) = Bz;
    CHxy(1:CPML,:,:) = bx(1:CPML,1).*CHxy(1:CPML,:,:) + axd(1:CPML,1).*diff(Ey(1:CPML+1,:,:),1,1);
    CHxy(CPML+1:end,:,:) = bx(CPML+1:end,1).*CHxy(CPML+1:end,:,:) + axd(CPML+1:end,1).*diff(Ey(Nx-CPML+1:Nx+1,:,:),1,1);
    CHyx(:,1:CPML,:) = by(1,1:CPML).*CHyx(:,1:CPML,:) + ayd(1,1:CPML).*diff(Ex(:,1:CPML+1,:),1,2);
    CHyx(:,CPML+1:end,:) = by(1,CPML+1:end).*CHyx(:,CPML+1:end,:) + ayd(1,CPML+1:end).*diff(Ex(:,Ny-CPML+1:Ny+1,:),1,2);
    Bz(:,:,:) = Bz(:,:,:) - kdx(1:Nx,1).*diff(Ey,1,1) + kdy(1,1:Ny).*diff(Ex,1,2);
    Bz(1:CPML,:,:) = Bz(1:CPML,:,:) - dt*CHxy(1:CPML,:,:);   Bz(Nx-CPML+1:Nx,:,:) = Bz(Nx-CPML+1:Nx,:,:) - dt*CHxy(CPML+1:end,:,:);
    Bz(:,1:CPML,:) = Bz(:,1:CPML,:) + dt*CHyx(:,1:CPML,:);   Bz(:,Ny-CPML+1:Ny,:) = Bz(:,Ny-CPML+1:Ny,:) + dt*CHyx(:,CPML+1:end,:);
    Hz = Hz + u_1.*( Bz - DT(1:Nx,1:Ny,:) );
    %% Electric Field Update
    Hsrc_x = Hx_srcArray(t)*G_Distrtn;
    Hsrc_y = Hy_srcArray(t)*G_Distrtn;
    il = ik-1;
    DT(1:Nx,:,:) = Dx;
    ET(1:Nx,:,:) = Ex;
    CDyz(:,1:CPML,:) = by(1,1:CPML).*CDyz(:,1:CPML,:) + ayd(1,1:CPML).*diff(Hz(:,1:CPML+1,:),1,2);
    CDyz(:,CPML+1:end,:) = by(1,CPML+1:end).*CDyz(:,CPML+1:end,:) + ayd(1,CPML+1:end).*diff(Hz(:,Ny-CPML:Ny,:),1,2);
    CDzy(:,:,1:CPML) = bz(1,1,1:CPML).*CDzy(:,:,1:CPML) + azd(1,1,1:CPML).*diff(Hy(:,:,1:CPML+1),1,3);
    CDzy(:,:,CPML+1:end) = bz(1,1,CPML+1:end).*CDzy(:,:,CPML+1:end) + azd(1,1,CPML+1:end).*diff(Hy(:,:,Nz-CPML:Nz),1,3);
    Dx(:,2:Ny,2:Nz) = Dx(:,2:Ny,2:Nz) + kdy(1,2:Ny).*diff(Hz(:,:,2:Nz),1,2) - kdz(1,1,2:Nz).*diff(Hy(:,2:Ny,:),1,3);
    %Dx(:,2:Ny,il) = Dx(:,2:Ny,il) + kdz(1,1,il).*Hsrc_y(1:Nx,2:Ny);  % TFSF
    Dx(:,1:CPML,:) = Dx(:,1:CPML,:) + dt*CDyz(:,1:CPML,:);   Dx(:,Ny-CPML+1:Ny,:) = Dx(:,Ny-CPML+1:Ny,:) + dt*CDyz(:,CPML+1:end,:);
    Dx(:,:,1:CPML) = Dx(:,:,1:CPML) - dt*CDzy(:,:,1:CPML);   Dx(:,:,Nz-CPML+1:Nz) = Dx(:,:,Nz-CPML+1:Nz) - dt*CDzy(:,:,CPML+1:end);
    Ex = s1(1:Nx,:,:).*Ex + s2(1:Nx,:,:).*( Dx - DT(1:Nx,:,:) ) + s3(1:Nx,:,:).*Jx;
    Jx = Kd(1:Nx,:,:).*Jx + Bd(1:Nx,:,:).*(Ex + ET(1:Nx,:,:));
    DT(:,1:Ny,:) = Dy;
    ET(:,1:Ny,:) = Ey;
    CDzx(:,:,1:CPML) = bz(1,1,1:CPML).*CDzx(:,:,1:CPML) + azd(1,1,1:CPML).*diff(Hx(:,:,1:CPML+1),1,3);
    CDzx(:,:,CPML+1:end) = bz(1,1,CPML+1:end).*CDzx(:,:,CPML+1:end) + azd(1,1,CPML+1:end).*diff(Hx(:,:,Nz-CPML:Nz),1,3);
    CDxz(1:CPML,:,:) = bx(1:CPML,1).*CDxz(1:CPML,:,:) + axd(1:CPML,1).*diff(Hz(1:CPML+1,:,:),1,1);
    CDxz(CPML+1:end,:,:) = bx(CPML+1:end,1).*CDxz(CPML+1:end,:,:) + axd(CPML+1:end,1).*diff(Hz(Nx-CPML:Nx,:,:),1,1);
    Dy(2:Nx,:,2:Nz) = Dy(2:Nx,:,2:Nz) + kdz(1,1,2:Nz).*diff(Hx(2:end-1,:,:),1,3) - kdx(2:Nx,1).*diff(Hz(:,:,2:Nz),1,1);
    %Dy(2:Nx,:,il) = Dy(2:Nx,:,il) - kdz(1,1,il).*Hsrc_x(2:Nx,1:Ny);  % TFSF
    Dy(:,:,1:CPML) = Dy(:,:,1:CPML) + dt*CDzx(:,:,1:CPML);   Dy(:,:,Nz-CPML+1:Nz) = Dy(:,:,Nz-CPML+1:Nz) + dt*CDzx(:,:,CPML+1:end);
    Dy(1:CPML,:,:) = Dy(1:CPML,:,:) - dt*CDxz(1:CPML,:,:);   Dy(Nx-CPML+1:Nx,:,:) = Dy(Nx-CPML+1:Nx,:,:) - dt*CDxz(CPML+1:end,:,:);
    Ey = s1(:,1:Ny,:).*Ey + s2(:,1:Ny,:).*( Dy - DT(:,1:Ny,:) ) + s3(:,1:Ny,:).*Jy;
    Jy = Kd(:,1:Ny,:).*Jy + Bd(:,1:Ny,:).*(Ey + ET(:,1:Ny,:));
    DT(:,:,1:Nz) = Dz;
    ET(:,:,1:Nz) = Ez;
    CDxy(1:CPML,:,:) = bx(1:CPML,1).*CDxy(1:CPML,:,:) + axd(1:CPML,1).*diff(Hy(1:CPML+1,:,:),1,1);
    CDxy(CPML+1:end,:,:) = bx(CPML+1:end,1).*CDxy(CPML+1:end,:,:) + axd(CPML+1:end,1).*diff(Hy(Nx-CPML:Nx,:,:),1,1);
    CDyx(:,1:CPML,:) = by(1,1:CPML).*CDyx(:,1:CPML,:) + ayd(1,1:CPML).*diff(Hx(:,1:CPML+1,:),1,2);
    CDyx(:,CPML+1:end,:) = by(1,CPML+1:end).*CDyx(:,CPML+1:end,:) + ayd(1,CPML+1:end).*diff(Hx(:,Ny-CPML:Ny,:),1,2);
    Dz(2:end-1,2:end-1,:) = Dz(2:end-1,2:end-1,:) + kdx(2:Nx,1).*diff(Hy(:,2:end-1,:),1,1) - kdy(1,2:Ny).*diff(Hx(2:end-1,:,:),1,2);
    Dz(1:CPML,:,:) = Dz(1:CPML,:,:) + dt*CDxy(1:CPML,:,:);   Dz(Nx-CPML+1:Nx,:,:) = Dz(Nx-CPML+1:Nx,:,:) + dt*CDxy(CPML+1:end,:,:);
    Dz(:,1:CPML,:) = Dz(:,1:CPML,:) - dt*CDyx(:,1:CPML,:);   Dz(:,Ny-CPML+1:Ny,:) = Dz(:,Ny-CPML+1:Ny,:) - dt*CDyx(:,CPML+1:end,:);
    Ez = s1(:,:,1:Nz).*Ez + s2(:,:,1:Nz).*( Dz - DT(:,:,1:Nz) ) + s3(:,:,1:Nz).*Jz;
    Jz = Kd(:,:,1:Nz).*Jz + Bd(:,:,1:Nz).*(Ez + ET(:,:,1:Nz));
    %% Point Source
    %Ez(round(Nx/2)+1,round(Ny/2)+1,round(Nz/2)) = Ex_srcArray(t);sin(2*pi*f0*dt*t)./(1+exp(-.2*(t-60)));
    %% Plotting
    if(mod(t,ntvis)==1)
        figure(f1);
        p1 = subplot(2,2,1);
        p1.XLabel.String = 'Ez field';
        slice(X,Y,Z,Ez,yslice,xslice,zslice);
        axis([0 y(end) 0 x(end) 0 z(end)]); xlabel Y; ylabel X; zlabel Z;
        FieldMax = 0.9*abs(max(max(max(abs(Ez)))))+1e-100;
        view([1 0.1 0.1]); clim([-FieldMax,FieldMax]);
        shading interp; colorbar; drawnow;
        p2 = subplot(2,2,2);
        p2.XLabel.String = 'Ex field';
        slice(X1,Y1,Z1,Ex,yslice,xslice,zslice);
        axis([0 y(end) 0 x(end) 0 z(end)]); xlabel Y; ylabel X; zlabel Z;
        FieldMax = 0.9*abs(max(max(max(abs(Ex)))))+1e-100;
        view([1 0.1 0.1]); clim([-FieldMax,FieldMax]);
        shading interp; colorbar; drawnow;
        p3 = subplot(2,2,3);
        p3.FontSize = 7;
        p3.XLabel.String = 'Ez field at y = 0.5 Ny';
        FieldMax = 0.5*abs(max(max(max(abs(Ez)))))+1e-100;
        imagesc(za,xa,squeeze(Ez(:,ceil(Ny/2),:)),[-FieldMax,FieldMax]);colorbar;
        p4 = subplot(2,2,4);
        p4.FontSize = 7;
        p4.XLabel.String = 'Ex field';
        FieldMax = 0.5*abs(max(max(max(abs(Ex)))))+1e-100;
        imagesc(za,xa,squeeze(Ex(:,ceil(Ny/2),:)),[-FieldMax,FieldMax]);colorbar;
        currentframe = getframe(f1);     % capture the frame to save movie
        writeVideo(v,currentframe);
    end
    waitbar((t/nt),WaitBar1,['Please wait... | ',num2str(t),'/',num2str(nt),' | ',num2str(floor((t/nt)*100)),'%']);
end
toc