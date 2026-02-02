clear
count=0;
QMIN=1;
H=.01;
A=2.0926e7;
GM=1.4077e16;
DISTKM=10000.;
ANGDEG=0.;
ANG=ANGDEG/57.3;
PHI=DISTKM*3280./A;
if QMIN==1;
    GAM=3.14159/4.-PHI/4.;
else
    GAM=25./57.3;
    GAM=20./57.3;
end
ALTKM=0.;
ALT=ALTKM*3280.;
R0=A+ALT;
TOP=GM*(1.-cos(PHI));
TEMP=R0*cos(GAM)/A-cos(PHI+GAM);
BOT=R0*cos(GAM)*TEMP;
V=sqrt(TOP/BOT);
VRX=V*cos(1.5708-GAM+ANG);
VRY=V*sin(1.5708-GAM+ANG);
XLAM=R0*V*V/GM;
TOP1=tan(GAM)*(1-cos(PHI))+(1-XLAM)*sin(PHI);
BOT1P=(1-cos(PHI))/(XLAM*cos(GAM)*cos(GAM));
BOT1=(2-XLAM)*(BOT1P+cos(GAM+PHI)/cos(GAM));
TOP2=2*cos(GAM);
BOT2=XLAM*((2/XLAM-1)^1.5);
TOP3=sqrt(2/XLAM-1);
BOT3=cos(GAM)/tan(PHI/2)-sin(GAM);
TEMP=(TOP2/BOT2)*atan2(TOP3,BOT3);
TF=R0*(TOP1/BOT1+TEMP)/(V*cos(GAM));
VIC=V/3280.;
SCOUNT=0.;
X=(A+ALT)*cos(ANG);
Y=(A+ALT)*sin(ANG);
XFIRST=X;
YFIRST=Y;
X1=VRX;
Y1=VRY;
T=0.;
while ~((T>10.) & ALT<0.)
     XOLD=X;
     YOLD=Y;
     X1OLD=X1;
     Y1OLD=Y1;
     STEP=1;
     FLAG=0;
     while STEP <=1
         if FLAG==1
             STEP=2;
             X=X+H*XD;
             Y=Y+H*YD;
             X1=X1+H*X1D;
             Y1=Y1+H*Y1D;
             T=T+H;
         end
         TEMBOT=(X^2+Y^2)^1.5;
         X1D=-GM*X/TEMBOT;
         Y1D=-GM*Y/TEMBOT;
         XD=X1;
         YD=Y1;
         ALT=sqrt(X^2+Y^2)-A;
         FLAG=1;
     end
     FLAG=0;
     X=(XOLD+X)/2+.5*H*XD;
     Y=(YOLD+Y)/2+.5*H*YD;
     X1=(X1OLD+X1)/2+.5*H*X1D;
     Y1=(Y1OLD+Y1)/2+.5*H*Y1D;
SCOUNT=SCOUNT+H;
if SCOUNT>=9.9999
     SCOUNT=0.;
     ALTKM=ALT/3280.;
     R=sqrt(X^2+Y^2);
     RF=sqrt(XFIRST^2+YFIRST^2);
     CBETA=(X*XFIRST+Y*YFIRST)/(R*RF);
     if CBETA>1
         CBETA=1;
     end
     BETA=acos(CBETA);
     DISTRTKM=A*BETA/3280.;
     VELKM=sqrt(X1^2+Y1^2)/3280.;
     count=count+1;
     ArrayT(count)=T;
     ArrayDISTRTKM(count)=DISTRTKM;
     ArrayALTKM(count)=ALTKM;
     ArrayVELKM(count)=VELKM;
end
end

figure
plot(ArrayDISTRTKM,ArrayALTKM),grid
xlabel('Downrange (km)')
ylabel('Altitude (km) ')
figure
plot(ArrayT,ArrayVELKM),grid
xlabel('Time (s)')
ylabel('Velocity (km/s) ')
clc
output=[ArrayT',ArrayDISTRTKM',ArrayALTKM',ArrayVELKM'];
save datfil.txt output -ascii
disp 'simulation finished'
