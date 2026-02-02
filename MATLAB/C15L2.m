clear
count=0;
H=.01
A=2.0926e7;
GM=1.4077e16;
GAM=45.;
ALTNM=0.;
V=24000.;
ANGDEG=0.;
ANG=ANGDEG/57.3;
VRX=V*cos(1.5708-GAM/57.3+ANG);
VRY=V*sin(1.5708-GAM/57.3+ANG);
ALT=ALTNM/6076.;
S=0.;
SCOUNT=0.;
R0=A+ALT;
R1=V*sin(GAM/57.3);
PSI=0.;
X=(A+ALT)*cos(ANG);
Y=(A+ALT)*sin(ANG);
XFIRST=X;
YFIRST=Y;
X1=VRX;
Y1=VRY;
T=0.;
while ALTNM > -.0001
 	R0OLD=R0;
	R1OLD=R1;
	PSIOLD=PSI;
	XOLD=X;
	YOLD=Y;
	X1OLD=X1;
	Y1OLD=Y1;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
			STEP=2;
 			R0=R0+H*R0D;
			R1=R1+H*R1D;
			PSI=PSI+H*PSID;
			X=X+H*XD;
			Y=Y+H*YD;
			X1=X1+H*X1D;
			Y1=Y1+H*Y1D;
			T=T+H;
		end
		PSID=(A+ALT)*V*cos(GAM/57.3)/(R0*R0);
		R1D=-GM/(R0*R0)+R0*PSID*PSID;
		R0D=R1;
		TEMBOT=(X^2+Y^2)^1.5;
		X1D=-GM*X/TEMBOT;
		Y1D=-GM*Y/TEMBOT;
		XD=X1;
		YD=Y1;
		FLAG=1;
	end
	FLAG=0;	
	R0=(R0OLD+R0)/2+.5*H*R0D;
	R1=(R1OLD+R1)/2+.5*H*R1D;
	PSI=(PSIOLD+PSI)/2+.5*H*PSID;
	X=(XOLD+X)/2+.5*H*XD;
	Y=(YOLD+Y)/2+.5*H*YD;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	Y1=(Y1OLD+Y1)/2+.5*H*Y1D;
	S=S+H;
	if S>=9.99999
		S=0.;
		SPOLARNM=A*PSI/6076.;
		ALTPOLARNM=(R0-A)/6076.;
		ALTNM=(sqrt(X^2+Y^2)-A)/6076.;
		R=sqrt(X^2+Y^2);
		RF=sqrt(XFIRST^2+YFIRST^2);
		CBETA=(X*XFIRST+Y*YFIRST)/(R*RF);
		BETA=acos(CBETA);
		DISTNM=A*BETA/6076.;
		count=count+1;
		ArrayT(count)=T;
		ArraySPOLARNM(count)=SPOLARNM;
		ArrayALTPOLARNM(count)=ALTPOLARNM;
		ArrayDISTNM(count)=DISTNM;
		ArrayALTNM(count)=ALTNM;
	end
end
figure
plot(ArraySPOLARNM,ArrayALTPOLARNM,ArrayDISTNM,ArrayALTNM),grid
xlabel('Downrange (Nmi)')
ylabel('Altitude (Nmi) ')
clc
output=[ArrayT',ArraySPOLARNM',ArrayALTPOLARNM',ArrayDISTNM',ArrayALTNM'];
save datfil.txt output -ascii
disp 'simulation finished'
	
