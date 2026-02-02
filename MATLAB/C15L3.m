clear
count=0;
H=.01;
A=2.0926e7;
GM=1.4077e16;
GAM=0.;
ALTNM=1000.;
ALT=ALTNM*6076.;
XLAM=1.;
V=sqrt(GM*XLAM/(A+ALT));
ANGDEG=90.;
ANG=ANGDEG/57.3;
VRX=V*cos(1.5708-GAM/57.3+ANG);
VRY=V*sin(1.5708-GAM/57.3+ANG);
S=0.;
SCOUNT=0.;
X=(A+ALT)*cos(ANG);
Y=(A+ALT)*sin(ANG);
XFIRST=X;
YFIRST=Y;
X1=VRX;
Y1=VRY;
T=0.;
TF=30000.;
while T < TF
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
	S=S+H;
	if S>=49.99999
		S=0.;
		XNM=X/6076.;
		YNM=Y/6076.;
		count=count+1;
		ArrayT(count)=T;
		ArrayXNM(count)=XNM;
		ArrayYNM(count)=YNM;
	end
end
figure
plot(ArrayXNM,ArrayYNM),grid
xlabel('X (Nmi)')
ylabel('Y (Nmi) ')
clc
output=[ArrayT',ArrayXNM',ArrayYNM'];
save datfil.txt output -ascii
disp 'simulation finished'
	
