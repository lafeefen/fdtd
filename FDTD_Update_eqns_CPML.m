%% Start loop
tic
for t = 1:nt
    %% Magnetic Field Update
    %Esrc_x = Ex_srcArray(t)*G_Distrtn;
    %Esrc_y = Ey_srcArray(t)*G_Distrtn;
    %ik = zSrc;
    DT(:,1:Ny,1:Nz) = Bx;
    CHyz = by(1,1:Ny).*CHyz + ayd(1,1:Ny).*diff(Ez,1,2);
    CHzy = bz(1,1,1:Nz).*CHzy + azd(1,1,1:Nz).*diff(Ey,1,3);
    Bx = Bx - kdy(1,1:Ny).*diff(Ez,1,2) + kdz(1,1,1:Nz).*diff(Ey,1,3) - dt*( CHyz - CHzy );
    %Bx(:,1:Ny,ik ) = Bx(:,1:Ny,ik) - kdz(1,1,ik).*Esrc_y(:,1:Ny);%  diff(Ey(:,1:Ny,zSrc:zSrc+1),1,3);  %TFSF
    Hx = Hx + u_1.*( Bx - DT(:,1:Ny,1:Nz) );
    DT(1:Nx,:,1:Nz)= By;
    CHzx = bz(1,1,1:Nz).*CHzx + azd(1,1,1:Nz).*diff(Ex,1,3);
    CHxz = bx(1:Nx,1).*CHxz + axd(1:Nx,1).*diff(Ez,1,1);
    By = By - kdz(1,1,1:Nz).*diff(Ex,1,3) + kdx(1:Nx,1).*diff(Ez,1,1) - dt*( CHzx - CHxz );
    %By(1:Nx,:,ik)  = By(1:Nx,:,ik) + kdz(1,1,ik).*Esrc_x(1:Nx,:);%  diff(Ex(:,:,zSrc:zSrc+1),1,3);     %TFSF
    Hy = Hy + u_1.*( By - DT(1:Nx,:,1:Nz) );
    DT(1:Nx,1:Ny,:) = Bz;
    CHxy = bx(1:Nx,1).*CHxy + axd(1:Nx,1).*diff(Ey,1,1);
    CHyx = by(1,1:Ny).*CHyx + ayd(1,1:Ny).*diff(Ex,1,2);
    Bz = Bz - kdx(1:Nx,1).*diff(Ey,1,1) + kdy(1,1:Ny).*diff(Ex,1,2) - dt*( CHxy - CHyx );
    Hz = Hz + u_1.*( Bz - DT(1:Nx,1:Ny,:) );
    %% Electric Field Update
    %Hsrc_x = Hx_srcArray(t)*G_Distrtn;
    %Hsrc_y = Hy_srcArray(t)*G_Distrtn;
    %il = ik-1;
    DT(1:Nx,:,:) = Dx;
    ET(1:Nx,:,:) = Ex;
    CDyz(:,2:Ny,2:Nz) = by(1,2:Ny).*CDyz(:,2:Ny,2:Nz) + ayd(1,2:Ny).*diff(Hz(:,:,2:Nz),1,2);
    CDzy(:,2:Ny,2:Nz) = bz(1,1,2:Nz).*CDzy(:,2:Ny,2:Nz) + azd(1,1,2:Nz).*diff(Hy(:,2:Ny,:),1,3);
    Dx(:,2:Ny,2:Nz) = Dx(:,2:Ny,2:Nz) + kdy(:,2:Ny).*diff(Hz(:,:,2:Nz),1,2) - kdz(:,:,2:Nz).*diff(Hy(:,2:Ny,:),1,3) + dt*( CDyz(:,2:Ny,2:Nz) - CDzy(:,2:Ny,2:Nz) );
    Ex = s1(1:Nx,:,:).*Ex + s2(1:Nx,:,:).*( Dx - DT(1:Nx,:,:) ) + s3(1:Nx,:,:).*Jx;
    Jx = Kd(1:Nx,:,:).*Jx + Bd(1:Nx,:,:).*(Ex + ET(1:Nx,:,:));
    DT(:,1:Ny,:) = Dy;
    ET(:,1:Ny,:) = Ey;
    CDzx(2:Nx,:,2:Nz) = bz(1,1,2:Nz).*CDzx(2:Nx,:,2:Nz) + azd(1,1,2:Nz).*diff(Hx(2:Nx,:,:),1,3);
    CDxz(2:Nx,:,2:Nz) = bx(2:Nx,1).*CDxz(2:Nx,:,2:Nz) + axd(2:Nx,1).*diff(Hz(:,:,2:end-1),1,1);
    Dy(2:Nx,:,2:Nz) = Dy(2:Nx,:,2:Nz) + kdz(1,1,2:Nz).*diff(Hx(2:Nx,:,:),1,3) - kdx(2:Nx,1).*diff(Hz(:,:,2:Nz),1,1) + dt*( CDzx(2:Nx,:,2:Nz) - CDxz(2:Nx,:,2:Nz) );
    Ey = s1(:,1:Ny,:).*Ey + s2(:,1:Ny,:).*( Dy - DT(:,1:Ny,:) ) + s3(:,1:Ny,:).*Jy;
    Jy = Kd(:,1:Ny,:).*Jy + Bd(:,1:Ny,:).*(Ey + ET(:,1:Ny,:));
    DT(:,:,1:Nz) = Dz;
    ET(:,:,1:Nz) = Ez;
    CDxy(2:Nx,2:Ny,:) = bx(2:Nx,1).*CDxy(2:Nx,2:Ny,:) + axd(2:Nx,1).*diff(Hy(:,2:end-1,:),1,1);
    CDyx(2:Nx,2:Ny,:) = by(1,2:Ny).*CDyx(2:Nx,2:Ny,:) + ayd(1,2:Ny).*diff(Hx(2:end-1,:,:),1,2);
    Dz(2:Nx,2:Ny,:) = Dz(2:Nx,2:Ny,:) + kdx(2:Nx,:).*diff(Hy(:,2:Ny,:),1,1) - kdy(:,2:Ny).*diff(Hx(2:Nx,:,:),1,2) + dt*( CDxy(2:Nx,2:Ny,:) - CDyx(2:Nx,2:Ny,:) );
    Ez = s1(:,:,1:Nz).*Ez + s2(:,:,1:Nz).*( Dz - DT(:,:,1:Nz) ) + s3(:,:,1:Nz).*Jz;
    Jz = Kd(:,:,1:Nz).*Jz + Bd(:,:,1:Nz).*(Ez + ET(:,:,1:Nz));
    %% Point Source
    Ez(round(Nx/2)+1,round(Ny/2)+1,round(Nz/2)) = sin(2*pi*f0*dt*t)./(1+exp(-.2*(t-60)));
    %% Plotting
    if(mod(t,ntvis)==0)
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
        waitbar((t/nt),WaitBar1,['Please wait... | ',num2str(t),'/',num2str(nt),' | ',num2str(floor((t/nt)*100)),'%']);
    end

end
toc