clear
count=0;
RDESKM=2000.;
% 1=IRBM,2=ICBM
ITGT=1;
TLOFT=200.;
TUPT=15.;
if ITGT==1
	TPZ=180.;
else
	TPZ=240.;
end
QMIN=1;
TF=2000.;
TFINISH=3000.;
LEFT=1;
QBOOST=1;
QOOMPH=1;
CW=0;
SWITCH=0;
GAMDEG=89.99;
H=.01;
T=0.;
S=0.;
A=2.0926E7;
GM=1.4077E16;
ALT=0;
ANGDEG=0.;
ANG=ANGDEG/57.3;
XLONGM=ANG;
X=(A+ALT)*cos(ANG);
Y=(A+ALT)*sin(ANG);
Z=0;
ALT=sqrt(X^2+Y^2)-A;
XFIRST=X;
YFIRST=Y;
ZFIRST=Z;
X1=cos(1.5708-GAMDEG/57.3+ANG);
Y1=sin(1.5708-GAMDEG/57.3+ANG);
AXT=0.;
AYT=0.;
XLONGTDEG=57.3*RDESKM*3280./A;
TF=252.+.223*RDESKM-(5.44E-6)*RDESKM*RDESKM;
TF=TF+TLOFT;
XLONGT=XLONGTDEG/57.3;
XF=A*cos(XLONGT);
YF=A*sin(XLONGT);
ZF=0;
while ALT>-1
	XOLD=X;
	YOLD=Y;
	X1OLD=X1;
	Y1OLD=Y1;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
        	STEP=2;
			X=X+H*X1;
			Y=Y+H*Y1;
			X1=X1+H*X1D;
			Y1=Y1+H*Y1D;
			T=T+H;
		end
		if ITGT==1
			if T<180.
				WGT=-212.*T+44000.;
				TRST=54100.;
			else
				WGT=3300.;
				TRST=0.;
			end
		else
			if T<120
				WGT=-2622*T+440660.;
				TRST=725850.;
			elseif T<240.
				WGT=-642.*T+168120.;
				TRST=182250.;
			else
				WGT=5500.;
				TRST=0.;
			end
		end
		AT=32.2*TRST/WGT;
		TEMBOT=(X^2+Y^2)^1.5;
		X1D=-GM*X/TEMBOT+AXT;
		Y1D=-GM*Y/TEMBOT+AYT;
		ALT=sqrt(X^2+Y^2)-A;
		FLAG=1;
	end
	FLAG=0;
	X=(XOLD+X)/2+.5*H*X1;
	Y=(YOLD+Y)/2+.5*H*Y1;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	Y1=(Y1OLD+Y1)/2+.5*H*Y1D;
 	S=S+H;
	if QBOOST==1
		TGOLAM=TF-T;
		[VRX,VRY,VRZ]=LAMBERT3D(X,Y,Z,TGOLAM,XF,YF,ZF,SWITCH);
		DELX=VRX-X1;
		DELY=VRY-Y1;
		DEL=sqrt(DELX^2+DELY^2);
		if T<TPZ & DEL>500
			AXT=AT*DELX/DEL;
			AYT=AT*DELY/DEL;
		elseif DEL<500
			TRST=0.;
			QBOOST=0;
			AXT=0.;
			AYT=0.;
			X1=VRX;
			Y1=VRY;
			X1OLD=X1;
			Y1OLD=Y1;
		else
			QBOOST=0;
			AXT=0.;
			AYT=0.;
		end
		if T<TUPT
			RTMAG=sqrt(X^2+Y^2);
			AXT=AT*X/RTMAG;
			AYT=AT*Y/RTMAG;
		end
	end
	if S>=.99999
 		S=0.;
 		DISTKM=distance3dkm(X,Y,Z,XFIRST,YFIRST,ZFIRST);
		ALTKM=(sqrt(X^2+Y^2)-A)/3280.;
		VELK=sqrt(X1^2+Y1^2)/3280.;
		count=count+1;
		ArrayT(count)=T;
		ArrayDISTKM(count)=DISTKM;
		ArrayALTKM(count)=ALTKM;
		ArrayVELK(count)=VELK;
	end
end
figure
plot(ArrayDISTKM',ArrayALTKM'),grid
xlabel('Downrange (km)')
ylabel('Altitude (km) ')
clc
output=[ArrayT',ArrayDISTKM',ArrayALTKM'];
save datfil.txt output -ascii
disp 'simulation finished'

	
