count=0;
PI=3.14159;
DEGRAD=360./(2.*PI);
H=.01;
A=2.0926e7;
GM=1.4077e16;
GAM=30.;
ALTNM=0.;
V=24000.;
TF=2200.;
ALT=ALTNM*6076.;
ANG=0.;
VRX=V*cos(PI/2.-GAM/DEGRAD+ANG);
VRY=V*sin(PI/2.-GAM/DEGRAD+ANG);
SCOUNT=0.;
X=(A+ALT)*cos(ANG);
Y=(A+ALT)*sin(ANG);
XFIRST=X;
YFIRST=Y;
X1=VRX;
Y1=VRY;
T=0.;
while T<(TF-.0001)
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
	X=(XOLD+X)/2+.5*H*XD;
	Y=(YOLD+Y)/2+.5*H*YD;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	Y1=(Y1OLD+Y1)/2+.5*H*Y1D;
	SCOUNT=SCOUNT+H;
	if SCOUNT>=1.99999
		SCOUNT=0.;
		ALTNM=(sqrt(X^2+Y^2)-A)/6076.;
		R=sqrt(X^2+Y^2);
		RF=sqrt(XFIRST^2+YFIRST^2);
		CBETA=(X*XFIRST+Y*YFIRST)/(R*RF);
		BETA=acos(CBETA);
		DISTNM=A*BETA/6076.;
		count=count+1;
		ArrayT(count)=T;
		ArrayDISTNM(count)=DISTNM;
		ArrayALTNM(count)=ALTNM;
	end
end
figure
plot(ArrayDISTNM,ArrayALTNM),grid
xlabel('Downrange (Nmi)')
ylabel('Altitude (Nmi) ')
clc
output=[ArrayT',ArrayDISTNM',ArrayALTNM'];
save datfil.txt output -ascii
X
Y
disp 'simulation finished'	