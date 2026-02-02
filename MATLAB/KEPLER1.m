function [X1]=KEPLER1(X0,T0,T1)
GMX=398923.;
REX=6380.;
TLIMIT=1.E-10;
KN=10;
MUQR = 1.;
DT = T1 - T0;
DX=10.;
if (abs(DT) <TLIMIT)
	for I=1:6
          X1(I) = X0(I);
    end

end
TIME_FACTOR = sqrt(REX^3/GMX);
VEX = REX/TIME_FACTOR;
DT = DT/TIME_FACTOR;
for I=1:3
         R0(I) = X0(I)/REX;
         V0(I) = X0(I+3)/VEX;
end
R0MAG = sqrt(R0(1)^2 +R0(2)^2 +R0(3)^2);
V0MAG = sqrt(V0(1)^2 +V0(2)^2 +V0(3)^2);
D0 = R0(1)*V0(1) +R0(2)*V0(2) +R0(3)*V0(3);
SIGMA0 = D0/MUQR;
ALP0 = 2./R0MAG - V0MAG*V0MAG;
if  ALP0 == 0.
        A0 = 1.E30;
else
        A0 = 1./ALP0;
end

X = ALP0*DT;
if  ALP0 <= 0. 
        X = .1*DT/R0MAG;
end
for K=1:KN
	if ALP0<0.
		Y = ALP0*X*X;
        	YQR = sqrt(-Y);
        	CY = (1. -cosh(YQR))/Y;
        	SY = (sinh(YQR) -YQR)/(YQR^3);
        elseif ALP0==0.
        	Y = 0.;
        	CY = .5;
        	SY = 1./6.;
        else
        	Y = ALP0*X*X;
        	YQR = sqrt(Y);;
        	CY = (1. -cos(YQR))/Y;
        	SY = (YQR -sin(YQR))/(YQR^3);
       end
     	U1 = X*(1.-Y*SY);
        U2 = X*X*CY;
        U3 = X*X*X*SY;
        FX = R0MAG*U1 +SIGMA0*U2 +U3 -DT*MUQR;
        DFX = SIGMA0*U1 +(1. -ALP0*R0MAG)*U2 +R0MAG;
        DFX2 = SIGMA0*(1. -Y*CY) +(1. -ALP0*R0MAG)*U1;
        SDFX = DFX/(abs(DFX));
        DX0 = 16.*DFX*DFX;
        DX1 = 20.*FX*DFX2;
        DX2 = 16.*DFX*DFX - 20.*FX*DFX2;
        if  DX2 > 0.
          DX = 5.*FX/(DFX +SDFX*sqrt(DX2));
        else
          DX = .5*X;
        end
         X =X -DX;
end

RMAG =DFX;
F = 1. -U2/R0MAG;
G = DT -U3/MUQR;
DF = -MUQR*U1/(RMAG*R0MAG);
DG = 1. -U2/RMAG;
for I=1:3
        X1(I) = (F*R0(I) +G*V0(I))*REX;
        X1(I+3) = (DF*R0(I) +DG*V0(I))*VEX;
end
