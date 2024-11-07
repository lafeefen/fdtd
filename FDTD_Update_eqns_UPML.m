tic
for t = 1:nt
    %% Magnetic Field Update
    GT(:,1:Ny,1:Nz) = Cx;
    DT(:,1:Ny,1:Nz) = Bx;
    Cx = sz1(1,1,1:Nz).*Cx + sz2(1,1,1:Nz).*( diff(Ey,1,3) - diff(Ez,1,2) );
    Bx = sy1(1,1:Ny).*Bx + sxy1(:,1:Ny).*Cx - sxy2(:,1:Ny).*GT(:,1:Ny,1:Nz);
    Hx = Hx + u_1.*( Bx - DT(:,1:Ny,1:Nz) );
    GT(1:Nx,:,1:Nz) = Cy;
    DT(1:Nx,:,1:Nz) = By;
    Cy = sx1(1:Nx,1).*Cy + sx2(1:Nx,1).*( diff(Ez,1,1) - diff(Ex,1,3) );
    By = sz1(1,1,1:Nz).*By + syz1(1,:,1:Nz).*Cy - syz2(1,:,1:Nz).*GT(1:Nx,:,1:Nz);
    Hy = Hy + u_1.*( By - DT(1:Nx,:,1:Nz) );
    GT(1:Nx,1:Ny,:) = Cz;
    DT(1:Nx,1:Ny,:) = Bz;
    Cz = sy1(1,1:Ny).*Cz + sy2(1,1:Ny).*( diff(Ex,1,2) - diff(Ey,1,1) );
    Bz = sx1(1:Nx,1).*Bz + szx1(1:Nx,1,:).*Cz - szx2(1:Nx,1,:).*GT(1:Nx,1:Ny,:);
    Hz = Hz + u_1.*( Bz - DT(1:Nx,1:Ny,:) );
    %% TFSF boundary
    % {
    Esrcx = Ex_srcArray(t);
    Esrcy = Ey_srcArray(t);
    ik = zSrc+1;
    for ii = 1:Nx
        for ij = 1:Ny
            %Ex(ii,ij,ik)
            Esrc_x = Esrcx*exp(-((((ii-xSrc).^2)+((ij-ySrc).^2))/xytau));    %exp(-(((Tstep-to)/tau).^2));   sin(2*pi*f0*dt*t)./(1+exp(-.2*(t-60)));
            Esrc_y = Esrcy*exp(-((((ii-xSrc).^2)+((ij-ySrc).^2))/xytau));    %exp(-(((Tstep-to)/tau).^2));   sin(2*pi*f0*dt*t)./(1+exp(-.2*(t-60)));
            Cx(ii,ij,ik) = Cx(ii,ij,ik) - sz2(1,1,ik)*( Esrc_y );
            Bx(ii,ij,ik) = Bx(ii,ij,ik) - sxy1(ii,ij)*sz2(1,1,ik)*( Esrc_y );
            Hx(ii,ij,ik) = Hx(ii,ij,ik) - u_1*sxy1(ii,ij)*sz2(1,1,ik)*( Esrc_y );
            Cy(ii,ij,ik) = Cy(ii,ij,ik) + sx2(ii,1)*( Esrc_x );
            By(ii,ij,ik) = By(ii,ij,ik) + syz1(1,ij,ik)*sx2(ii,1)*( Esrc_x );
            Hy(ii,ij,ik) = Hy(ii,ij,ik) + u_1*syz1(1,ij,ik)*sx2(ii,1)*( Esrc_x );
        end
    end
    %}
    %% Electric Field Update
    GT(1:Nx,:,:) = Gx;
    DT(1:Nx,:,:) = Dx;
    Gx(:,2:end-1,2:end-1) = sz1(1,1,2:end-1).*Gx(:,2:end-1,2:end-1) + sz2(1,1,2:end-1).*( diff(Hz(:,:,2:end-1),1,2) - diff(Hy(:,2:end-1,:),1,3) );
    Dx = sy1.*Dx + sxy1(1:Nx,:).*Gx - sxy2(1:Nx,:).*GT(1:Nx,:,:);
    GT(1:Nx,:,:) = Ex;
    Ex = s1(1:Nx,:,:).*Ex + s2(1:Nx,:,:).*( Dx - DT(1:Nx,:,:) ) + s3(1:Nx,:,:).*Jx;
    Jx = Kd(1:Nx,:,:).*Jx + Bd(1:Nx,:,:).*(Ex + GT(1:Nx,:,:));
    GT(:,1:Ny,:) = Gy;
    DT(:,1:Ny,:) = Dy;
    Gy(2:end-1,:,2:end-1) = sx1(2:end-1,1).*Gy(2:end-1,:,2:end-1) + sx2(2:end-1,1).*( diff(Hx(2:end-1,:,:),1,3) - diff(Hz(:,:,2:end-1),1,1) );
    Dy = sz1.*Dy + syz1(1,1:Ny,:).*Gy - syz2(1,1:Ny,:).*GT(:,1:Ny,:);
    GT(:,1:Ny,:) = Ey;
    Ey = s1(:,1:Ny,:).*Ey + s2(:,1:Ny,:).*( Dy - DT(:,1:Ny,:) ) + s3(:,1:Ny,:).*Jy;
    Jy = Kd(:,1:Ny,:).*Jy + Bd(:,1:Ny,:).*(Ey + GT(:,1:Ny,:));
    GT(:,:,1:Nz) = Gz;
    DT(:,:,1:Nz) = Dz;
    Gz(2:end-1,2:end-1,:) = sy1(1,2:end-1).*Gz(2:end-1,2:end-1,:) + sy2(1,2:end-1).*( diff(Hy(:,2:end-1,:),1,1) - diff(Hx(2:end-1,:,:),1,2) );
    Dz = sx1.*Dz + szx1(:,1,1:Nz).*Gz - szx2(:,1,1:Nz).*GT(:,:,1:Nz);
    GT(:,:,1:Nz) = Ez;
    Ez = s1(:,:,1:Nz).*Ez + s2(:,:,1:Nz).*( Dz - DT(:,:,1:Nz) ) + s3(:,:,1:Nz).*Jz;
    Jz = Kd(:,:,1:Nz).*Jz + Bd(:,:,1:Nz).*(Ez + GT(:,:,1:Nz));
    %% TFSF boundary
    % {
    Hsrcx = Hx_srcArray(t);
    Hsrcy = Hy_srcArray(t);
    ik = zSrc-1;
    for ii = 1:Nx
        for ij = 1:Ny
            Hsrc_x = Hsrcx*exp(-((((ii-xSrc).^2)+((ij-ySrc).^2))/xytau));    %exp(-(((Tstep-to)/tau).^2));   sin(2*pi*f0*dt*t)./(1+exp(-.2*(t-60)));
            Hsrc_y = Hsrcy*exp(-((((ii-xSrc).^2)+((ij-ySrc).^2))/xytau));    %exp(-(((Tstep-to)/tau).^2));   sin(2*pi*f0*dt*t)./(1+exp(-.2*(t-60)));
            Gx(ii,ij,ik) = Gx(ii,ij,ik) - sz2(1,1,ik)*( Hsrc_y );
            Dx(ii,ij,ik) = Dx(ii,ij,ik) - sxy1(ii,ij)*sz2(1,1,ik)*( Hsrc_y );
            Ex(ii,ij,ik) = Ex(ii,ij,ik) - s2(ii,ij,ik)*sxy1(ii,ij)*sz2(1,1,ik)*( Hsrc_y );
            Jx(ii,ij,ik) = Jx(ii,ij,ik) - Bd(ii,ij,ik)*s2(ii,ij,ik)*sxy1(ii,ij)*sz2(1,1,ik)*( Hsrc_y );
            Gy(ii,ij,ik) = Gy(ii,ij,ik) + sx2(ii,1)*( Hsrc_x );
            Dy(ii,ij,ik) = Dy(ii,ij,ik) + syz1(1,ij,ik)*sx2(ii,1)*( Hsrc_x );
            Ey(ii,ij,ik) = Ey(ii,ij,ik) + s2(ii,ij,ik)*syz1(1,ij,ik)*sx2(ii,1)*( Hsrc_x );
            Jy(ii,ij,ik) = Jy(ii,ij,ik) + Bd(ii,ij,ik)*s2(ii,ij,ik)*syz1(1,ij,ik)*sx2(ii,1)*( Hsrc_x );
        end
    end
    %}
    %% Point Source
    %{
    Esrc_x = Ex_srcArray(t);
    ik = round(Nz/2)+1;
    xSrc = round(Nx/2); ySrc = round(Ny/2);
    for ii = 1:Nx
        for ij = 1:Ny
            Ex(ii,ij,ik) = Esrc_x*exp(-((((ii-xSrc).^2)+((ij-ySrc).^2))/200));    %exp(-(((Tstep-to)/tau).^2));   sin(2*pi*f0*dt*t)./(1+exp(-.2*(t-60)));
        end
    end
    %}
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
        FieldMax = 0.9*abs(max(max(max(abs(Ez)))))+1e-100;
        imagesc(za,xa,squeeze(Ez(:,ceil(Ny/2),:)),[-FieldMax,FieldMax]);colorbar;
        p4 = subplot(2,2,4);
        p4.FontSize = 7;
        p4.XLabel.String = 'Ex field';
        FieldMax = 0.9*abs(max(max(max(abs(Ex)))))+1e-100;
        imagesc(za,xa,squeeze(Ex(:,ceil(Ny/2),:)),[-FieldMax,FieldMax]);colorbar;
    end
end
toc