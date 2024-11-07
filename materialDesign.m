function [InMat] = materialDesign(InMat,Nx,Ny,Nz,CPML)%,dx)
    rD0 = 40;             % radius of sphere in units of dx;
    cDx = ceil(Nx/2);     % position of Dielectric sphere's center in x-direction
    cDy = ceil(Ny/2);     % position of Dielectric sphere's center in x-direction
    cDz = CPML+100;       % position of Dielectric sphere's center in x-direction
    %{
    for ii = 1:Nx
        for ij = 1:Ny
            for ik = 1:Nz
                %if ((ii-cDx)^2 + (ij-cDy)^2 + (ik-cDz)^2 <= (rD0)^2) && ik <= (cDz-ceil(rD0/2))% && (L(1)-ii*dx<=ij*dy))
                if (((ii-cDx)^2 + (ij-cDy)^2) > (rD0-6)^2 ) && (((ii-cDx)^2 + (ij-cDy)^2) < (rD0)^2 ) && ik <= (cDz-ceil(rD0/2)+6)
                    InMat(ii,ij,ik) = 5;
                end
            end
        end
    end
    %}
    %{
    for ii = 1:Nx
        for ij = 1:Ny
            for ik = 1:Nz
                if ((ii-cDx)^2 + (ij-cDy)^2 + (ik-cDz)^2 <= (rD0)^2) && ik <= (cDz-ceil(rD0/2))% && (L(1)-ii*dx<=ij*dy))
                    InMat(ii,ij,ik) = 2;
                end
            end
        end
    end
    %}
    cDz = CPML+60;       % position of Dielectric sphere's center in x-direction
    %{
    for ii = 1:Nx
        for ij = 1:Ny
            for ik = 1:Nz
                if ((ii-cDx)^2 + (ij-cDy)^2 + (ik-cDz)^2 <= (rD0)^2) && ik >= (cDz+ceil(rD0/2))% && (L(1)-ii*dx<=ij*dy))
                    InMat(ii,ij,ik) = 2;
                end
            end
        end
    end
    %}
    cDx = ceil(Nx/2);     % position of Dielectric sphere's center in x-direction
    cDy = ceil(Ny/2);     % position of Dielectric sphere's center in x-direction
    cDz = CPML+160; %ceil(Nz/2)-50;  % position of Dielectric sphere's center in x-direction
    %{
    for ii = 1:Nx
        for ij = 1:Ny
            for ik = 1:Nz
                if ((ii-cDx)^2 + (ij-cDy)^2 + (ik-cDz)^2 <= (rD0)^2) && ik <= (cDz-ceil(rD0/2))% && (L(1)-ii*dx<=ij*dy))
                    InMat(ii,ij,ik) = 2;
                end
            end
        end
    end
    %}
    %{
    for ii = 1:Nx
        for ij = 1:Ny
            for ik = (cDz-ceil(rD0/2)+1):(cDz-ceil(rD0/2)+120)
                InMat(ii,ij,ik) = 2;
            end
        end
    end
    %}
    %{
    for ii = 1:Nx
        for ij = 1:Ny
            for ik = (cDz-ceil(rD0/2)+121):(cDz-ceil(rD0/2)+125)
                InMat(ii,ij,ik) = 3;
            end
        end
    end
    %}
    imagesc(squeeze(InMat(:,Ny/2,:)));
    %InMat((CPML+5):(Nx-CPML-4),(CPML+5):(Ny-CPML-4),CPML+10:CPML+20)=3;
    %InMat((CPML+5):(Nx-CPML-4),(CPML+5):(Ny-CPML-4),CPML+21:CPML+50)=2;
end