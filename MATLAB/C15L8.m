clear
n=0;
QFIRST=1;
H=.01;
GM=1.4077E16;
A=2.0926E7;
GAMDEG=170.;
GAM=GAMDEG/57.3;
DISTKM=10000.;
ANGDEG=0. ;
ANG=ANGDEG/57.3;
PHI=DISTKM*3280./A;
ALTKM=0.;
ALT=ALTKM*3280.;
R0=A+ALT;
TOP=GM*(1.-cos(PHI));
TEMP=R0*cos(GAM)/A-cos(PHI+GAM);
BOT=R0*cos(GAM)*TEMP;
V=sqrt(TOP/BOT);
VRX=V*cos(1.5708-GAM+ANG);
VRY=V*sin(1.5708-GAM+ANG);
S=0.;
SCOUNT=0.;
X=(A+ALT)*cos(ANG);
Y=(A+ALT)*sin(ANG);
XFIRST=X;
YFIRST=Y;
X1=VRX;
Y1=VRY;
IPZ=1;
T=0.;
BETAOLD=0.;
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
		FLAG=1;
     	end
     	FLAG=0;
	Y=(YOLD+Y)/2+.5*H*YD;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	Y1=(Y1OLD+Y1)/2+.5*H*Y1D;
	ALT=sqrt(X^2+Y^2)-A; 
	S=S+H;
	SCOUNT=SCOUNT+H;
	if SCOUNT>=9.9999
     		SCOUNT=0.;
		ALTKM=ALT/3280.;
		VKM=sqrt(X1^2+Y1^2)/3280.;
		R=sqrt(X^2+Y^2);
		RF=sqrt(XFIRST^2+YFIRST^2);
		CBETA=(X*XFIRST+Y*YFIRST)/(R*RF);
		if CBETA>1
			CBETA=1;
		end
		BETA=acos(CBETA);
		if IPZ==0
			BETA=2*3.14159-BETA;
		end
		DISTRTKM=A*BETA/3280;
		if T>100 & BETA<BETAOLD & QFIRST==1
			IPZ=0;
			QFIRST=0
		elseif QFIRST==1
			IPZ=1 ;
		end
		XKM=X/3280.; 
		YKM=Y/3280.;
		BETAOLD=BETA;
		n=n+1;
		ArrayT(n)=T;
		ArrayDISTRTKM(n)=DISTRTKM;
		ArrayALTKM(n)=ALTKM;
		ArrayXKM(n)=XKM;
		ArrayYKM(n)=YKM;
		ArrayVKM(n)=VKM;
	end
end
figure
plot(ArrayDISTRTKM,ArrayALTKM),grid
xlabel('Downrange (km)')
ylabel('Altitude (km)')
clc
output=[ArrayDISTRTKM',ArrayALTKM',ArrayVKM'];
save datfil.txt output -ascii
disp 'simulation finished'

