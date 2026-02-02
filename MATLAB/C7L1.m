clear
count=0;
VC=4000;
XNT=96.6;
YIC=0;
VM=3000;
HEDEG=0;
BETA=0.3;
XNP=3;
SIGNOISE=.001;
TF=10;
TS=.1;
NOISE=1;
Y=YIC;
YD=-VM*HEDEG/57.3;
YDIC=YD;
T=0;
H=.01;
S=0;
GFILTER=1.-BETA^2;
HFILTER=(1.-BETA)^2;
XLAMH=0;
XLAMDH=0;
XNC=0;
while T <= (TF - 1e-5)
 	YOLD=Y;
	YDOLD=YD;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
       			Y=Y+H*YD;
 			YD=YD+H*YDD;
			T=T+H;
			STEP=2;
		end;
     		TGO=TF-T+1e-5;
		RTM=VC*TGO;
		XLAM=Y/(VC*TGO);
		XLAMD=(RTM*YD+Y*VC)/(RTM^2);
		YDD=XNT-XNC;
		FLAG=1;
	end; 
	FLAG=0;
	Y=.5*(YOLD+Y+H*YD);
 	YD=.5*(YDOLD+YD+H*YDD);
	S=S+H;
 	if S>=(TS - 1e-5)
		S=0.;
		if NOISE==1
			XLAMNOISE=SIGNOISE*randn;
		else
			XLAMNOISE=0.;
		end;
		RES=XLAM-(XLAMH+TS*XLAMDH)+XLAMNOISE;
		XLAMH=GFILTER*RES+XLAMH+TS*XLAMDH;
		XLAMDH=HFILTER*RES/TS+XLAMDH;
		XNC=XNP*VC*XLAMDH;
		count=count+1;
		ArrayT(count)=T;
		ArrayY(count)=Y;
		ArrayXNC(count)=XNC;
		ArrayXLAMD(count)=XLAMD;
		ArrayXLAMDH(count)=XLAMDH;
 	end; 
end; 
figure
plot(ArrayT,ArrayXLAMD,ArrayT,ArrayXLAMDH),grid
title('Decreasing beta increase noise transmission of fading memory filter')
xlabel('Time (S)')
ylabel('Line of Sight Rate (Rad/S) ')
axis([0 10 -.01 .06])
output=[ArrayT',ArrayY',ArrayXNC',ArrayXLAMD',ArrayXLAMDH'];
save datfil.txt output -ascii
disp('Simulation Complete')






