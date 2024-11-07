function [Ex_srcArray,Ey_srcArray,Hx_srcArray,Hy_srcArray] = FDTD_Source(Source, STEPS, dt, Lamda, e0, u0, c0,Px)
    fundamentalConstantsUnits;
    T = 1:STEPS;
    Tstep = T*dt;
    delt = 1.5*dt;
    f0 = c/Lamda;        % center frequency of source excitation (Hz)
    Omega = 2*pi*f0;     % center frequency in radians
    %Ez(round(Nx/2)+1,round(Ny/2)+1,round(Nz/2)) = sin(2*pi*f0*dt*t)./(1+exp(-.2*(t-60)));
    Py = sqrt(1-(Px^2));
    if Source == 1
        tau = 4*fs; %2*fs;
        to = 8*fs;
        A =  sqrt(1);
        Ex_srcArray = Px*exp(-(((Tstep-to)/tau).^2));
        Ey_srcArray = Py*exp(-(((Tstep-to)/tau).^2));
        Hx_srcArray = -Py*A*exp(-(((Tstep-to + delt)/tau).^2))/(c*Mu0); %plot(T,Esrc,T,Hsrc);
        Hy_srcArray = Px*A*exp(-(((Tstep-to + delt)/tau).^2))/(c*Mu0); %plot(T,Esrc,T,Hsrc);
        if(STEPS>7000)
            Ex_srcArray(7000:STEPS) = 0;
            Ey_srcArray(7000:STEPS) = 0;
            Hx_srcArray(7000:STEPS) = 0;
            Hy_srcArray(7000:STEPS) = 0;
        end
    end
    if Source == 2
        tau = 1*fs; %3D: 3*fs; %2*fs; 1D:
        to = 3*tau; %3D: 10*fs; 1D:
        A = sqrt(1)/c0; %3D: sqrt(1); 1D: 
        Ex_srcArray = Px*(exp(-(((Tstep-to)/tau).^2)/2).*sin(2*pi*f0*(Tstep - to)));%/(c*Ep0);
        Ey_srcArray = Py*(exp(-(((Tstep-to)/tau).^2)/2).*sin(2*pi*f0*(Tstep - to)));
        Hx_srcArray = -Py*A*exp(-(((Tstep-to + delt)/tau).^2)/2).*sin(2*pi*f0*(Tstep - to + delt));%/(c*Mu0);%plot(T,Esrc,T,Hsrc);
        Hy_srcArray = Px*A*exp(-(((Tstep-to + delt)/tau).^2)/2).*sin(2*pi*f0*(Tstep - to + delt));
        if(STEPS>600)
            Ex_srcArray(600:STEPS) = 0;
            Ey_srcArray(600:STEPS) = 0;
            Hx_srcArray(600:STEPS) = 0;
            Hy_srcArray(600:STEPS) = 0;
        end
    end
    if Source == 3      %for Ricker Wavelet
        tor = 2*fs;
        A = sqrt(1)/c0; %3D: sqrt(1); 1D: 
        Ex_srcArray = Px*(1-2*(pi*f0*(Tstep-tor)).^2).*exp(-(pi*f0*(Tstep-tor)).^2);
        Ey_srcArray = Py*(1-2*(pi*f0*(Tstep-tor)).^2).*exp(-(pi*f0*(Tstep-tor)).^2);
        Hx_srcArray = -Py*A*((1-2*(pi*f0*(Tstep-tor+delt)).^2).*exp(-(pi*f0*(Tstep-tor+delt)).^2))/(c*Mu0);
        Hy_srcArray = Px*A*((1-2*(pi*f0*(Tstep-tor+delt)).^2).*exp(-(pi*f0*(Tstep-tor+delt)).^2))/(c*Mu0);
        if(STEPS>7000)
            Ex_srcArray(7000:STEPS) = 0;
            Ey_srcArray(7000:STEPS) = 0;
            Hx_srcArray(7000:STEPS) = 0;
            Hy_srcArray(7000:STEPS) = 0;
        end
    end
end
