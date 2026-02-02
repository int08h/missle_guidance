clear
count=0;
SWITCH1=0;
XLONGTDEG=7.42;
XLATTDEG=43.75;
XLATFDEG=36.175;
XLONGFDEG=-115.136;
TF=2000.;
A=2.0926E7;
GM=1.4077E16;
W=-6.283185/86400.;
XLONGF=XLONGFDEG/57.3;
XLATF=XLATFDEG/57.3;
XLONGT=XLONGTDEG/57.3;
XLATT=XLATTDEG/57.3;
XLONGF=XLONGF-W*TF;
XF=A*cos(XLATF)*cos(XLONGF);
YF=A*cos(XLATF)*sin(XLONGF);
ZF=A*sin(XLATF);
XT=A*cos(XLATT)*cos(XLONGT);
YT=A*cos(XLATT)*sin(XLONGT);
ZT=A*sin(XLATT);
[VRX,VRY,VRZ]=LAMBERT3D(XT,YT,ZT,TF,XF,YF,ZF,SWITCH1);
XTD=VRX;
YTD=VRY;
ZTD=VRZ;
XTINIT=XT;
YTINIT=YT;
ZTINIT=ZT;
T=0.;
H=.001;
S=0.;
ALTNM=(sqrt(XT^2+YT^2+ZT^2)-A)/6076.;
while ALTNM>-1
 	XTOLD=XT;
	YTOLD=YT;
	ZTOLD=ZT;
	XTDOLD=XTD;
	YTDOLD=YTD;
	ZTDOLD=ZTD;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
         		STEP=2;
			XT=XT+H*XTD;
			YT=YT+H*YTD;
			ZT=ZT+H*ZTD;
			XTD=XTD+H*XTDD;
			YTD=YTD+H*YTDD;
			ZTD=ZTD+H*ZTDD;
        		T=T+H;
        	end
        	TEMPBOTT=(XT^2+YT^2+ZT^2)^1.5;
		XTDD=-GM*XT/TEMPBOTT;
		YTDD=-GM*YT/TEMPBOTT;
		ZTDD=-GM*ZT/TEMPBOTT;
		ALTNM=(sqrt(XT^2+YT^2+ZT^2)-A)/6076.;
		FLAG=1;
	end
	FLAG=0;
	XT=.5*(XTOLD+XT+H*XTD);
 	YT=.5*(YTOLD+YT+H*YTD);
	ZT=.5*(ZTOLD+ZT+H*ZTD);
	XTD=.5*(XTDOLD+XTD+H*XTDD);
 	YTD=.5*(YTDOLD+YTD+H*YTDD);
	ZTD=.5*(ZTDOLD+ZTD+H*ZTDD);
	S=S+H;
	if S>=9.9999
		S=0.;
		XTE=XT*cos(W*T)-YT*sin(W*T);
		YTE=XT*sin(W*T)+YT*cos(W*T);
		ZTE=ZT;
		XLATT= atan2(ZTE, sqrt(XTE^2+YTE^2));
		XLATTDEG=57.3*XLATT;
		XLONGT= atan2(YTE, XTE);
		XLONGTDEG=57.3*XLONGT;
		DISTRTNM=distance3d(XTE,YTE,ZTE,XTINIT,YTINIT,ZTINIT);
		count=count+1;
		ArrayT(count)=T;
		ArrayDISTRTNM(count)=DISTRTNM;
		ArrayALTNM(count)=ALTNM;
		ArrayXLONGTDEG(count)=XLONGTDEG;
		ArrayXLATTDEG(count)=XLATTDEG;
	end
end

figure
plot(ArrayDISTRTNM,ArrayALTNM),grid
xlabel('Downrange (Nmi)')
ylabel('Altitude (Nmi) ')
clc
output=[ArrayT',ArrayDISTRTNM',ArrayALTNM'];
save datfil.txt output -ascii
output=[ArrayT',ArrayXLONGTDEG',ArrayXLATTDEG'];
csvwrite('trajfil.txt',output)

disp 'simulation finished'
	
