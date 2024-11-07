%% ===> FDTD with Drude Model <===
clear; close all; clc;tic
%% === Define constants and parameters <===
THz = 10^12; PHz = 10^15;
e0 = 8.854*10^-12;                     % Permittivity of vacuum   [F/m]
u0 = 4*pi*10^-7;                       % Permeability of vacuum   [H/m]
c0 = 1/(e0*u0)^.5;                     % Speed of light           [m/s]
Lamda = 785*10^-9;                     % Wavelength               [nm]
f0    = c0/Lamda;%1e6;                 % Frequency of Source      [Hz]
L0    = c0/f0;                         % Freespace Wavelength     [m]
t0    = 1/f0;                          % Source Period            [s]
Lf         = 100;%157;                 % Divisions per Wavelength
[Lx,Ly,Lz] = deal(1,1,1);              % Dimension of the system in wavelengths x,y,z (4.8,4.8,7.0);
nt         = 2000;                     % Number of time steps
ntvis      = 10;                       % time step for visualization
Nx = round(Lx*Lf);                     % Number of Nodes in x
Ny = round(Ly*Lf);                     % Number of Nodes in y
Nz = round(Lz*Lf);                     % Number of Nodes in z
x  = linspace(0,Lx,Nx+1)*L0;           % x vector                    [m]
y  = linspace(0,Ly,Ny+1)*L0;           % y vector                    [m]
z  = linspace(0,Lz,Nz+1)*L0;           % z vector                    [m]
[dx,dy,dz] = deal(x(2),y(2),z(2));     % x,y,z increment             [m]
dt = (dx^-2+dy^-2+dz^-2)^-.5/c0*0.99;  % Time step for CFL condition [s]
CPML = 5;                              % Number of PML nodes in the grid

% Drude model parameters %if Model == 7 % Drude Model CPML   %Au 9.1eV 13.8PHz 0.011PHz   %Ag 9.2eV 14.0PHz 0.032PHz ref:
WpAuD = 13.80*PHz;                           % Gold - Plasmon frequency    [Hz]
GammaD_Au = 11.0*THz;                        % Gold - Relaxation rate      [Hz]
WpAgD = 14.0*PHz;                            % Silver - Plasmon frequency  [Hz]   13.280*PHz;
GammaD_Ag = 32.0*THz;                        % Silver - Relaxation rate    [Hz]   91.269*THz;
wD = [0 0 WpAuD WpAgD 0];                    % Plasmon frequency           [Hz]
gD = [0 0 GammaD_Au GammaD_Ag 0];            % Relaxation rate             [Hz]
NoOfMedia = 5;                               % No. of different media: Vacuum, Silica, Gold, Silver, conductor
EpInf_Au  = 1.4447;                          % Gold
EpInf_Ag  = 1.4447;                          % Silver
eR = [1.0 1.5111 EpInf_Au EpInf_Ag 1.0];     % Permittivity of media [Vac, SiO2, Au, Ag]
s0 = [0.0 0.0 63.2*10^6 63.2*10^6 100*10^9]; % Conductivity of media
MurMedia = [1.0 1.0 1.0 1.0 1.0];            % Permeability of media
%% ===> Initialize Magnetic and Electric Field Vectors and Coefficients <===
Pml = 2; % 1: UPML, 2: CPML
Sep = 1; % 1: same, 2: Seperate - fdtd update equations for pml
if Pml == 1
    [GT,Gx,Gy,Gz,DT,Dx,Dy,Dz,Jx,Jy,Jz,Ex,Ey,Ez,Cx,Cy,Cz,Bx,By,Bz,Hx,Hy,Hz,G_Distrtn] = intialize_E_H_Fields_UPML(Nx,Ny,Nz);
elseif Pml == 2
    if Sep == 1
        [DT,Dx,Dy,Dz,Jx,Jy,Jz,ET,Ex,Ey,Ez,Bx,By,Bz,Hx,Hy,Hz,CDxy,CDyx,CDyz,CDzy,CDzx,CDxz,CHxy,CHyx,CHyz,CHzy,CHzx,CHxz,G_Distrtn] = intialize_E_H_Fields_CPML(Nx,Ny,Nz);
    elseif Sep == 2
        [DT,Dx,Dy,Dz,Jx,Jy,Jz,ET,Ex,Ey,Ez,Bx,By,Bz,Hx,Hy,Hz,CDxy,CDyx,CDyz,CDzy,CDzx,CDxz,CHxy,CHyx,CHyz,CHzy,CHzx,CHxz,bx,by,bz,axd,ayd,azd,kdx,kdy,kdz,G_Distrtn] = intialize_E_H_Fields_CPML_sep(Nx,Ny,Nz,CPML);
    end
end
%% ===> Design material <===
InMat = deal(ones(Nx+1,Ny+1,Nz+1));  % material space
[InMat] = materialDesign(InMat,Nx,Ny,Nz,CPML);
%% ===> H & E Coeffcients <===
[Kd,Bd,s1,s2,s3,u_1] = materialCoefficient(e0,u0,eR,gD,wD,dt,s0,InMat);
%% ===> PML Coeffcients <===
if Pml == 1
    [sx1,sy1,sz1,sx2,sy2,sz2,sxy1,sxy2,syz1,syz2,szx1,szx2] = UPML_Region_Sep(CPML,dx,dy,dz,dt,Nx,Ny,Nz,e0,u0);
elseif Pml == 2
    [kdx,kdy,kdz,bx,by,bz,axd,ayd,azd,bx1,by1,bz1,axd1,ayd1,azd1] = CPML_Region(3,Sep,CPML,dx,dy,dz,dt,Nx,Ny,Nz);
end
%%  ===> FDTD Source <===
Source = 2;     %1: Gaussian envelope 2: Modulated Gaussian 3: Ricker Wavelet
xSrc = round(Nx/2); ySrc = round(Ny/2); zSrc = CPML+ 10; % Position of the source
xytau=100; % width of spatial Gaussian envelope
[Ex_srcArray,Ey_srcArray,Hx_srcArray,Hy_srcArray] = FDTD_Source(Source, nt, dt, Lamda, e0, u0, c0,1);
for ii = 1:Nx
    for ij = 1:Ny
        G_Distrtn(ii,ij) = exp(-((((ii-xSrc).^2)+((ij-ySrc).^2))/xytau));
    end
end
%%  ===> Start loop <===
xa = linspace(1,Nx,Nx);
ya = linspace(1,Ny+1,Ny+1);
za = linspace(1,Nz+1,Nz+1);
f1 = figure('Position', [70 30 1400 700]);%f1 = figure(1);%f1 = figure('Renderer', 'painters', 'Position', [10 10 1400 400]);
set(gcf, 'Renderer', 'zbuffer');
p1 = subplot(2,2,1);
[X,Y,Z] = meshgrid(y,x,z(1:end-1));     % for plotting Ez
xslice = x(end)/2;                      % location of y-z planes
yslice = y(end)/2;                      % location of x-z planes
zslice = z(end)/2;                      % location of x-y planes
p2 = subplot(2,2,2);
[X1,Y1,Z1] = meshgrid(y,x(1:end-1),z);  % for plotting Ez
p3 = subplot(2,2,3);
p4 = subplot(2,2,4);
v = VideoWriter(['FDTD_v6_Nx_',num2str(Nx),'_CPML_',num2str(CPML),'_ADE.avi']);
open(v);
WaitBar1 = waitbar(0,['Please wait... | ','0%']);
toc
if Pml == 1
    FDTD_Update_eqns_UPML;
elseif Pml == 2
    if Sep == 1
        FDTD_Update_eqns_CPML;
    elseif Sep == 2
        FDTD_Update_eqns_CPML_sep;
    end
end
close(WaitBar1);
close(v);
disp("Done!!!");
sound(sin(1:3000));