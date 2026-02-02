clear
n=0;
XNTG=0.;
HEDEG=0.;
XNP=3.;
RM1IC=0.;
RM2IC=10000.;
RT1IC=30000.;
RT2IC=0.;
VM=3000.;
VT=0.;
XNCLIMG=9999999.;
APN=1;
XLAMFDEG=-90.;
H=.0001;
XNCLIM=32.2*XNCLIMG;
XLAMF=XLAMFDEG/57.3;
XNT=32.2*XNTG;
RM1=RM1IC;
RM2=RM2IC;
RT1=RT1IC;
RT2=RT2IC;
BETA=0.;
VT1=-VT*cos(BETA);
VT2=VT*sin(BETA);
HE=HEDEG/57.3;
T=0.;
S=0.;
RTM1=RT1-RM1;
RTM2=RT2-RM2;
RTM=sqrt(RTM1^2+RTM2^2);
XLAM=atan2(RTM2,RTM1); 
XLEAD=asin(VT*sin(BETA+XLAM)/VM);
THET=XLAM+XLEAD;
VM1=VM*cos(THET+HE);
VM2=VM*sin(THET+HE);
VTM1=VT1-VM1;
VTM2=VT2-VM2;
VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
while VC >= 0 
 	if RTM < 1000
		H=.00001;
	else
		H=.0001;
	end
 	BETAOLD=BETA;
	RT1OLD=RT1;
	RT2OLD=RT2;
	RM1OLD=RM1;
	RM2OLD=RM2;
	VM1OLD=VM1;
	VM2OLD=VM2;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
         		STEP=2;
 			BETA=BETA+H*BETAD;
 			RT1=RT1+H*VT1;
			RT2=RT2+H*VT2;
			RM1=RM1+H*VM1;
			RM2=RM2+H*VM2;
			VM1=VM1+H*AM1;
			VM2=VM2+H*AM2;
			T=T+H;
		end
		RTM1=RT1-RM1;
		RTM2=RT2-RM2;
		RTM=sqrt(RTM1^2+RTM2^2);
		VTM1=VT1-VM1;
		VTM2=VT2-VM2;
		VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
		XLAM=atan2(RTM2,RTM1);
		XLAMD=(RTM1*VTM2-RTM2*VTM1)/(RTM*RTM);
		TGO=RTM/VC;
		if APN==0
			XNC=XNP*VC*XLAMD;
		else
			XNT1=XNT*sin(BETA);
			XNT2=XNT*cos(BETA);
			XNTPLOS=-XNT1*sin(XLAM)+XNT2*cos(XLAM);
			XNC=4.*VC*XLAMD+XNTPLOS+2.*VC*(XLAM-XLAMF)/TGO;
		end
		if XNC>XNCLIM
			XNC=XNCLIM;
		end
		if XNC<-XNCLIM
			XNC=-XNCLIM;
		end
		AM1=-XNC*sin(XLAM);
		AM2=XNC*cos(XLAM);
		VT1=-VT*cos(BETA);
		VT2=VT*sin(BETA);
		if VT==0.
			BETAD=0.;
		else
			BETAD=XNT/VT;
		end
		FLAG=1;
	end
	FLAG=0;
 	BETA=.5*(BETAOLD+BETA+H*BETAD);
 	RT1=.5*(RT1OLD+RT1+H*VT1);
	RT2=.5*(RT2OLD+RT2+H*VT2);
	RM1=.5*(RM1OLD+RM1+H*VM1);
	RM2=.5*(RM2OLD+RM2+H*VM2);
	VM1=.5*(VM1OLD+VM1+H*AM1);
	VM2=.5*(VM2OLD+VM2+H*AM2);
	S=S+H;
	if S>=.09999	
		S=0.;
		n=n+1;
		RT1K=RT1/1000.;
		RT2K=RT2/1000.;
		RM1K=RM1/1000.;
		RM2K=RM2/1000.;
		XLAMDEG=XLAM*57.3;
		XNCG=XNC/32.2;
		ArrayT(n)=T;
		ArrayRT1K(n)=RT1K;
		ArrayRT2K(n)=RT2K;
		ArrayRM1K(n)=RM1K;
		ArrayRM2K(n)=RM2K;
		ArrayXNCG(n)=XNCG;
		ArrayXLAMDEG(n)=XLAMDEG;
	end
end
RTM
figure
plot(ArrayRT1K,ArrayRT2K,ArrayRM1K,ArrayRM2K),grid
title('Engagement Geometry')
xlabel('Downrange (Kft) ')
ylabel('Altitude (Kft)')
figure
plot(ArrayT,ArrayXNCG),grid
title('Commanded Acceleration')
xlabel('Time (Sec) ')
ylabel('XNC (G)')
axis([0 14 -20 25])
figure
plot(ArrayT,ArrayXLAMDEG),grid
title('Line-of-Sight Angle')
xlabel('Time (Sec) ')
ylabel('XLAM (Deg)')
axis([0 14 -100 0])
clc
output=[ArrayT',ArrayRT1K',ArrayRT2K',ArrayRM1K',...
		ArrayRM2K',ArrayXNCG',ArrayXLAMDEG'];
save datfil.txt output -ascii
disp '*** Simulation Complete'
	
