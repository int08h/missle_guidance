clear
count=0;
VC=4000.;
XNT=96.6;
YIC=0.;
VM=3000.;
HEDEG=0.;
BETA=.8;
XNP=3.;
SIGNOISE=.001;
TF=10.;
TS=.1;
NOISE=1;
Y=YIC;
YD=-VM*HEDEG/57.3;
YDIC=YD;
T=0.;
H=.01;
S=0.;
GFILTER=1.-BETA^3;
HFILTER=1.5*((1.-BETA)^2)*(1.+BETA);
KFILTER=.5*((1.-BETA)^3);
YH=0.;
YDH=0.;
XNTH=0.;
XNC=0.;
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
      		end
      		TGO=TF-T+.00001;
      		RTM=VC*TGO;
      		XLAM=Y/(VC*TGO);
      		XLAMD=(RTM*YD+Y*VC)/(RTM^2);
      		YDD=XNT-XNC;
      		FLAG=1;
	end
	FLAG=0;
 	Y=.5*(YOLD+Y+H*YD);
 	YD=.5*(YDOLD+YD+H*YDD);
	S=S+H;
	if S>=(TS - 1e-5)
      		S=0.;
      		if NOISE==1,
			XLAMNOISE=SIGNOISE*randn;
		else
			XLAMNOISE=0.;
      		end;	
      		YSTAR=RTM*(XLAM+XLAMNOISE);
      		RES=YSTAR-YH-TS*YDH-.5*TS*TS*(XNTH-XNC);
      		YH=GFILTER*RES+YH+TS*YDH+.5*TS*TS*(XNTH-XNC);
      		YDH=HFILTER*RES/TS+YDH+TS*(XNTH-XNC);
      		XNTH=2.*KFILTER*RES/(TS*TS)+XNTH;
      		XLAMDH=(YH+YDH*TGO)/(VC*TGO*TGO);
      		XNC=XNP*VC*XLAMDH;
      		count=count+1;
      		ArrayT(count)=T;
      		ArrayY(count)=Y;
      		ArrayXNCG(count)=XNC/32.2;
      		ArrayXLAMD(count)=XLAMD;
      		ArrayXLAMDH(count)=XLAMDH;
      		ArrayXNTG(count)=XNT/32.2;
      		ArrayXNTHG(count)=XNTH/32.2;
   	end
end
figure
plot(ArrayT,ArrayXLAMD,ArrayT,ArrayXLAMDH),grid
xlabel('Time (S)')
ylabel('Line of Sight Rate (Rad/S) ')
axis([0 10 0 .05])
figure
plot(ArrayT,ArrayXNTG,ArrayT,ArrayXNTHG),grid
xlabel('Time (S)')
ylabel('Acceleration (G) ')
clc
output=[ArrayT',ArrayY',ArrayXNCG',ArrayXLAMD',ArrayXLAMDH',ArrayXNTG',ArrayXNTHG'];
save datfil.txt output -ascii
disp 'simulation finished'



