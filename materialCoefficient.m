function [Kd,Bd,s1,s2,s3,u_1] = materialCoefficient(e0,u0,eR,gD,wD,dt,s0,InMat)
    Kd = ( (2 - gD(InMat).*dt)./(2 + gD(InMat).*dt) );
    Bd = ( (e0*dt.*wD(InMat).^2)./(2 + gD(InMat).*dt) );
    s1 = ( (2*e0.*eR(InMat) - s0(InMat).*dt - Bd(InMat).*dt)./(2*e0.*eR(InMat) + s0(InMat).*dt + Bd(InMat).*dt) );
    s2 = ( 2./(2*e0.*eR(InMat) + s0(InMat).*dt + Bd(InMat).*dt) );
    s3 = ( ((1+Kd(InMat)).*dt)./(2*e0.*eR(InMat) + s0(InMat).*dt + Bd(InMat).*dt) );
    u_1 = 1/u0;
end