clear
count=0;
RDESKM=7000.;
TF=2000.;
TFINISH=999999.;
TLOFT=500.;
TGRAVEND=100.;
GAMDEGIC=89.8;
TUPT=20.;
RDESRKM=560.;
LEFT=1;
QBOOST=1;
QOOMPH=1;
QFINISH=1;
QZERO=0;
CW=0;	
SWITCH=0;
QFIRST=1;
GAMDEG=GAMDEGIC;
H=.01;
T=0.;
S=0.;
A=2.0926E7;
GM=1.4077E16;
ALTNM=0.;
ALT=ALTNM*6076.;
ANGDEG=0.;
ANG=ANGDEG/57.3;
XLONGM=ANG;
X=(A+ALT)*cos(ANG);
Y=(A+ALT)*sin(ANG);
ALT=sqrt(X^2+Y^2)-A;
XFIRST=X;
YFIRST=Y;
X1=cos(1.5708-GAMDEG/57.3+ANG);
Y1=sin(1.5708-GAMDEG/57.3+ANG);
AXT=0.;
AYT=0.;
XLONGTDEG=57.3*RDESKM*3280./A;
XLONGRDEG=57.3*RDESRKM*3280./A;
TF=252.+.223*RDESKM-(5.44E-6)*RDESKM*RDESKM;
TF=TF+TLOFT;
XLONGT=XLONGTDEG/57.3;
XLONGR=XLONGRDEG/57.3;
XF=A*cos(XLONGT);
YF=A*sin(XLONGT);
XR=A*cos(XLONGR);
YR=A*sin(XLONGR);
AXT=0.;
AYT=0.;
Z=0;
ZF=0;
ZFIRST=0;
while ~(ALTNM<-1 | T>TFINISH)
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
		end;
		if T<120.
			WGT=-2622*T+440660.;
			TRST=725850.;
		elseif T<240.
			WGT=-642.*T+168120.;
			TRST=182250.;
		else
			WGT=5500.;
			TRST=0.;
		end
		AT=32.2*TRST/WGT;
		XD=X1;
		YD=Y1;
		VEL=sqrt(XD^2+YD^2);
		TEMBOT=(X^2+Y^2)^1.5;
		X1D=-GM*X/TEMBOT+AXT;
		Y1D=-GM*Y/TEMBOT+AYT;
		ALT=sqrt(X^2+Y^2)-A;
		ACCG=sqrt(AXT^2+AYT^2)/32.2;
		FLAG=1;
	end
	FLAG=0;
	X=(XOLD+X)/2+.5*H*XD;
	Y=(YOLD+Y)/2+.5*H*YD;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	Y1=(Y1OLD+Y1)/2+.5*H*Y1D;
	S=S+H;
	Z=0.;
	ZF=0.;
	TGOLAM=TF-T;
	if (QBOOST==1 & T>TGRAVEND)
		TGOLAM=TF-T;
		[VRX,VRY,VRZ]=LAMBERT3D(X,Y,Z,TGOLAM,XF,YF,ZF,SWITCH);
		DELX=VRX-X1;
		DELY=VRY-Y1;
		DEL=sqrt(DELX^2+DELY^2);
		if (T<240 & DEL>500.)
			AXT=AT*DELX/DEL;
			AYT=AT*DELY/DEL;
		elseif DEL<500.
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
	elseif (T>=TUPT & T<= TGRAVEND & QFIRST==1)
		QFIRST=0;
		VEL=sqrt(XD^2+YD^2);
		X1=VEL*cos(1.5708-GAMDEGIC/57.3+ANG);
		Y1=VEL*sin(1.5708-GAMDEGIC/57.3+ANG);
		X1OLD=X1;
		Y1OLD=Y1;
		AXT=AT*X1/VEL;
		AYT=AT*Y1/VEL;
	elseif (T>=TUPT & T<=TGRAVEND)
		VEL=sqrt(XD^2+YD^2);
		AXT=AT*X1/VEL;
		AYT=AT*Y1/VEL;
	elseif T<=TUPT
		RTMAG=sqrt(X^2+Y^2);
		AXT=AT*X/RTMAG;
		AYT=AT*Y/RTMAG;
	end
	if S>=.9999
		S=0.;
		DISTNM=distance3dkm(X,Y,Z,XFIRST,YFIRST,ZFIRST);
		ALTNM=(sqrt(X^2+Y^2)-A)/3280.;
		RMAG=sqrt(X^2+Y^2);
		VMAG=sqrt(XD^2+YD^2);
		GAMDEG=90-57.3*acos((X*XD+Y*YD)/(RMAG*VMAG));
		RHO=.0034*exp(-ALT/22000.);
		Q=.5*RHO*VEL*VEL;
		RRMAG=sqrt(XR^2+YR^2);
		RRTMAG=sqrt((X-XR)^2+(Y-YR)^2);
		ELDEG=90.-57.3*acos((XR*(X-XR)+YR*(Y-YR))/(RRMAG*RRTMAG));
		if (ELDEG>2. & ELDEG<85)
			ISEE=1;
		else
			ISEE=0;
		end
		count=count+1;
		ArrayT(count)=T;
		ArrayDISTNM(count)=DISTNM;
		ArrayALTNM(count)=ALTNM;
		ArrayQ(count)=Q;
		ArrayISEE(count)=ISEE;
	end
end
figure
plot(ArrayDISTNM,ArrayALTNM),grid
xlabel('Downrange (km)')
ylabel('Altitude (km) ')
clc
output=[ArrayT',ArrayDISTNM',ArrayALTNM',ArrayQ',ArrayISEE'];
save datfil.txt output -ascii
disp 'simulation finished'
