clear
count=0;
H=.01;	
A=2.0926e7;
GM=1.4077e16;
GAMDEG=23.;
GAM=GAMDEG/57.3;
DISTNM=6000.;
ANGDEG=0.;
ANG=ANGDEG/57.3;
PHI=DISTNM*6076./A;
ALTNM=0.;
ALT=ALTNM*6076.;
R0=A+ALT;
TOP=GM*(1.-cos(PHI));
TEMP=R0*cos(GAM)/A-cos(PHI+GAM);
BOT=R0*cos(GAM)*TEMP;
V=sqrt(TOP/BOT);
VRX=V*cos(1.5708-GAM+ANG);
VRY=V*sin(1.5708-GAM+ANG);
S=0.;
X=(A+ALT)*cos(ANG);
Y=(A+ALT)*sin(ANG);
XFIRST=X;
YFIRST=Y;
X1=VRX;
Y1=VRY;
T=0.;
while ALT >=0.
	XOLD=X;
	YOLD=Y;
	X1OLD=X1;
	Y1OLD=Y1;
	STEP=1;
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
	X=(XOLD+X)/2+.5*H*XD;
	Y=(YOLD+Y)/2+.5*H*YD;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	Y1=(Y1OLD+Y1)/2+.5*H*Y1D;
	S=S+H;
	if S>=9.99999
		S=0.;
		XNM=X/6076.;
		YNM=Y/6076.;
		ALT=sqrt(X^2+Y^2)-A;
		ALTNM=ALT/6076.;
		R=sqrt(X^2+Y^2);
		RF=sqrt(XFIRST^2+YFIRST^2);
		CBETA=(X*XFIRST+Y*YFIRST)/(R*RF);
		BETA=acos(CBETA);
		DISTNM=A*BETA/6076.;
		count=count+1;
		ArrayT(count)=T;
		ArrayXNM(count)=XNM;
		ArrayYNM(count)=YNM;
		ArrayDISTNM(count)=DISTNM;
		ArrayALTNM(count)=ALTNM;
	end
end
figure
plot(ArrayDISTNM,ArrayALTNM),grid
xlabel('Downrange (Nmi)')
ylabel('Altitude (Nmi) ')
clc
output=[ArrayT',ArrayXNM',ArrayYNM',ArrayDISTNM',ArrayALTNM'];
save datfil.txt output -ascii
disp 'simulation finished'
	
