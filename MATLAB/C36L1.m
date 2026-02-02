clear
n=0;
XNT=0.;
HEDEG=-20.;
XNCLIM=999999.;
PN=0;
XLAMFDEG=-30.;
VC=4000.;
VM=3000.;
TF=10.;
XNP=3.;
XLAMF=XLAMFDEG/57.3;
Y=0.;
YD=-VM*HEDEG/57.3;
T=0.;
H=.001;
S=0.;
while T<=(TF-.0001)
 	YOLD=Y;
	YDOLD=YD;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
         		STEP=2;
 			Y=Y+H*YD;
 			YD=YD+H*YDD;
			T=T+H;
		end
		TGO=TF-T+.00001;
 		XLAM=Y/(VC*TGO);
		XLAMD=(Y+YD*TGO)/(VC*TGO*TGO);
		if PN==1
			XNC=XNP*VC*XLAMD;
		else
			XNC=4.*VC*XLAMD+XNT+2.*VC*(XLAM-XLAMF)/TGO;
		end
		if XNC>XNCLIM
			XNC=XNCLIM;
		end
		if XNC<-XNCLIM
			XNC=-XNCLIM;
		end
		YDD=XNT-XNC;
		FLAG=1;
	end
	FLAG=0;
 	Y=.5*(YOLD+Y+H*YD);
 	YD=.5*(YDOLD+YD+H*YDD);
	S=S+H;
	if S>=.09999	
		S=0.;
		n=n+1;
		XLAMDEG=XLAM*57.3;
		XNCG=XNC/32.2;
		ArrayT(n)=T;
		ArrayY(n)=Y;
		ArrayXNCG(n)=XNCG;
		ArrayXLAMDEG(n)=XLAMDEG;
	end
end
figure
plot(ArrayT,ArrayY),grid
title('Relative Trajectory')
xlabel('Time (Sec) ')
ylabel('Y (Ft)')
figure
plot(ArrayT,ArrayXNCG),grid
title('Commanded Acceleration')
xlabel('Time (Sec) ')
ylabel('XNC (G)')
axis([0 10 -40 30])
figure
plot(ArrayT,ArrayXLAMDEG),grid
title('Line-of-Sight Angle')
xlabel('Time (Sec) ')
ylabel('XLAM (Deg)')
axis([0 10 -30 10])
clc
output=[ArrayT',ArrayY',ArrayXNCG',ArrayXLAMDEG'];
save datfil.txt output -ascii
disp '*** Simulation Complete'
