clear all
close all
count=0;
TS=.001;
XK3=-.362;
TA=.831;
B11=.00461;
B12=.00136;
XKFB=.00134;
WZFB=395.;
Z1=.015;
W1=259.;
XKFB2=.000664;
WZFB2=255.;
Z2=.022;
W2=649.;
WT=W1;
ZT=Z1;
WT2=W2;
ZT2=Z2;
WB=W1;
ZB=.5;
WB2=W2;
ZB2=.5;
TF=1.;
XIN=1.;
WACT=100.;
WACT=400.
ZACT=.7;
XKR=.1
XKR=.3
% Notch=1 (first mode filter), Notch=2 (both filters),Notch=3 (no filters)
NOTCH=3;
NOTCH=1;
NOTCH=2;
XK=(1.-XKR*XK3)/(XKR*XK3)
% FB1=1. (First mode) FB2=1 (Second mode)
FB1=1.;
FB2=1.;
T=0.;
S=0.;
H=.00001;
E=0.;
ED=0.;
E1=0.;
E1D=0.;
E1DD=0.;
DEL=0.;
DELD=0.;
DELC=0;
E2=0.;
E2D=0.;
E3=0.;
E3D=0.;
E3DD=0.;
E4=0.;
E4D=0.;
THDFIL1=0;
THDFIL=0;
while T<=TF
    EOLD=E;
    EDOLD=ED;
    E1OLD=E1;
    E1DOLD=E1D;
    E1DDOLD=E1DD;
    DELOLD=DEL;
    DELDOLD=DELD;
    E2OLD=E2;
    E2DOLD=E2D;
    E3OLD=E3;
    E3DOLD=E3D;
    E3DDOLD=E3DD;
    E4OLD=E4;
    E4DOLD=E4D;

    DELDD=WACT*WACT*(DELC-DEL-2.*ZACT*DELD/WACT);
    EDD=(DEL-E-B11*ED)/B12;
    E1DDD=W1*W1*(DELD-E1D-2.*Z1*E1DD/W1);
    E3DDD=W2*W2*(DELD-E3D-2.*Z2*E3DD/W2);
    THDRB=XK3*(E+TA*ED);
    THDFB=XKFB*(E1D+E1DDD/WZFB^2);
    THDFB2=XKFB2*(E3D+E3DDD/WZFB2^2);
    THDTOT=THDRB+FB1*THDFB+FB2*THDFB2;
    if NOTCH==1
        DELC=XKR*(XK*XIN+THDFIL1);
    elseif NOTCH==2
        DELC=XKR*(XK*XIN+THDFIL);
    else
        DELC=XKR*(XK*XIN+THDTOT);
    end
    E2DD=WB*WB*(THDTOT-E2-2*ZB*E2D/WB);
    THDFIL1=E2+2*ZT*E2D/WT+E2DD/WT^2;
    E4DD=WB2*WB2*(THDFIL1-E4-2*ZB2*E4D/WB2);
    THDFIL=E4+2*ZT2*E4D/WT2+E4DD/WT2^2;
    
    E=E+H*ED;
    ED=ED+H*EDD;
    E1=E1+H*E1D;
    E1D=E1D+H*E1DD;
    E1DD=E1DD+H*E1DDD;
    DEL=DEL+H*DELD;
    DELD=DELD+H*DELDD;
    E2=E2+H*E2D;
    E2D=E2D+H*E2DD;
    E3=E3+H*E3D;
    E3D=E3D+H*E3DD;
    E3DD=E3DD+H*E3DDD;
    E4=E4+H*E4D;
    E4D=E4D+H*E4DD;
    T=T+H;

    DELDD=WACT*WACT*(DELC-DEL-2.*ZACT*DELD/WACT);
    EDD=(DEL-E-B11*ED)/B12;
    E1DDD=W1*W1*(DELD-E1D-2.*Z1*E1DD/W1);
    E3DDD=W2*W2*(DELD-E3D-2.*Z2*E3DD/W2);
    THDRB=XK3*(E+TA*ED);
    THDFB=XKFB*(E1D+E1DDD/WZFB^2);
    THDFB2=XKFB2*(E3D+E3DDD/WZFB2^2);
    THDTOT=THDRB+FB1*THDFB+FB2*THDFB2;
    if NOTCH==1
        DELC=XKR*(XK*XIN+THDFIL1);
    elseif NOTCH==2
        DELC=XKR*(XK*XIN+THDFIL);
    else
        DELC=XKR*(XK*XIN+THDTOT);
    end
    E2DD=WB*WB*(THDTOT-E2-2*ZB*E2D/WB);
    THDFIL1=E2+2*ZT*E2D/WT+E2DD/WT^2;
    E4DD=WB2*WB2*(THDFIL1-E4-2*ZB2*E4D/WB2);
    THDFIL=E4+2*ZT2*E4D/WT2+E4DD/WT2^2;

    E=.5*(EOLD+E+H*ED);
    ED=.5*(EDOLD+ED+H*EDD);
    E1=.5*(E1OLD+E1+H*E1D);
    E1D=.5*(E1DOLD+E1D+H*E1DD);
    E1DD=.5*(E1DDOLD+E1DD+H*E1DDD);
    DEL=.5*(DELOLD+DEL+H*DELD);
    DELD=.5*(DELDOLD+DELD+H*DELDD);
    E2=.5*(E2OLD+E2+H*E2D);
    E2D=.5*(E2DOLD+E2D+H*E2DD);
    E3=.5*(E3OLD+E3+H*E3D);
    E3D=.5*(E3DOLD+E3D+H*E3DD);
    E3DD=.5*(E3DDOLD+E3DD+H*E3DDD);
    E4=.5*(E4OLD+E4+H*E4D);
    E4D=.5*(E4DOLD+E4D+H*E4DD);

    S=S+H;
    if S>=(TS-.00001)
        S=0.;
        count=count+1;
        ArrayT(count)=T;
        ArrayTHDTOT(count)=THDTOT;
	end
end
figure
plot(ArrayT,ArrayTHDTOT),grid
xlabel('Time (Sec)')
ylabel('Total Body Rate (deg/s)')
clc
output=[ArrayT',ArrayTHDTOT'];
save datfil.txt output  -ascii
disp 'simulation finished'
