clear all
close all
count=0;
TS=.001
XK3=-.362;
TA=.831;
B11=.00461;
B12=.00136;
XKFB=.00134;
WZFB=395.;
Z1=.015;
W1=259.;
TF=1.;
DELC=1.;
WACT=100.;
ZACT=.7;
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
while T<=TF
    EOLD=E;
    EDOLD=ED;
    E1OLD=E1;
    E1DOLD=E1D;
    E1DDOLD=E1DD;
    DELOLD=DEL;
    DELDOLD=DELD;
    DELDD=WACT*WACT*(DELC-DEL-2.*ZACT*DELD/WACT);
    EDD=(DEL-E-B11*ED)/B12;
    E1DDD=W1*W1*(DELD-E1D-2.*Z1*E1DD/W1);
    THDRB=XK3*(E+TA*ED);
    THDFB=XKFB*(E1D+E1DDD/WZFB^2);
    THDTOT=THDRB+THDFB;
    E=E+H*ED;
    ED=ED+H*EDD;
    E1=E1+H*E1D;
    E1D=E1D+H*E1DD;
    E1DD=E1DD+H*E1DDD;
    DEL=DEL+H*DELD;
    DELD=DELD+H*DELDD;
    T=T+H;
    DELDD=WACT*WACT*(DELC-DEL-2.*ZACT*DELD/WACT);
    EDD=(DEL-E-B11*ED)/B12;
    E1DDD=W1*W1*(DELD-E1D-2.*Z1*E1DD/W1);
    THDRB=XK3*(E+TA*ED);
    THDFB=XKFB*(E1D+E1DDD/WZFB^2);
    THDTOT=THDRB+THDFB;
    E=.5*(EOLD+E+H*ED);
    ED=.5*(EDOLD+ED+H*EDD);
    E1=.5*(E1OLD+E1+H*E1D);
    E1D=.5*(E1DOLD+E1D+H*E1DD);
    E1DD=.5*(E1DDOLD+E1DD+H*E1DDD);
    DEL=.5*(DELOLD+DEL+H*DELD);
    DELD=.5*(DELDOLD+DELD+H*DELDD);
    S=S+H;
    if S>=(TS-.00001)
        S=0.;
        count=count+1;
        ArrayT(count)=T;
        ArrayTHDRB(count)=THDRB;
        ArrayTHDFB(count)=THDFB;
        ArrayTHDTOT(count)=THDTOT;
    end
end
figure
plot(ArrayT,ArrayTHDRB),grid
xlabel('Time (Sec)')
ylabel('Rigid Body Rate (deg/s)')
figure
plot(ArrayT,ArrayTHDFB),grid
xlabel('Time (Sec)')
ylabel('Flexible Body Rate (deg/s)')
figure
plot(ArrayT,ArrayTHDTOT),grid
xlabel('Time (Sec)')
ylabel('Total Body Rate (deg/s)')
clc
output=[ArrayT',ArrayTHDRB',ArrayTHDFB',ArrayTHDTOT'];
save datfil.txt output  -ascii
disp 'simulation finished'
