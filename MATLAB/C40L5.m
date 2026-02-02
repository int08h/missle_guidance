clear
n=0;
XNTG=4.;
VT=1000.;
VM=3000.;
RM1=0.;
RM2=10000.;
RM3=-1000.;
RT1=30000.;
RT2=10000.;
RT3=0.;
GAMFPDEG= -30.;
GAMFYDEG= 20.;
XNT=32.2*XNTG;
BETA=0.;
VT1=-VT*cos(BETA);
VT2=VT*sin(BETA);
VT3=0;
GAMFP=GAMFPDEG/57.3;
GAMFY=GAMFYDEG/57.3;
[VM1,VM2,VM3,TF]=LAUNCHLOGIC(RM1,RM2,RM3,RT1,RT2,RT3,VT1,VT2,...
							 VT3,VM);
H=.0001;
T=0.;
S=0.;
RTM1=RT1-RM1;
RTM2=RT2-RM2;
RTM3=RT3-RM3;
RTM=sqrt(RTM1^2+RTM2^2+RTM3^2);
VTM1=VT1-VM1;
VTM2=VT2-VM2;
VTM3=VT3-VM3;
VC=-(RTM1*VTM1+RTM2*VTM2+RTM3*VTM3)/RTM;
VM1F=VM*cos(GAMFP)*cos(GAMFY);
VM2F=VM*sin(GAMFP);
VM3F=VM*cos(GAMFP)*sin(GAMFY);
while ~(VC<0)
 	if RTM<1000.
		H=.00001;
	else
		H=.0001; 
	end
	BETAOLD=BETA;
	RT1OLD=RT1;
	RT2OLD=RT2;
	RT3OLD=RT3;
	RM1OLD=RM1;
	RM2OLD=RM2;
	RM3OLD=RM3;
	VM1OLD=VM1;
	VM2OLD=VM2;
	VM3OLD=VM3;
	STEP=1;
	FLAG=0;
	while STEP<=1
		if FLAG==1
			STEP=2;
			BETA=BETA+H*BETAD;
			RT1=RT1+H*VT1;
			RT2=RT2+H*VT2;
			RT3=RT3+H*VT3;
			RM1=RM1+H*VM1;
			RM2=RM2+H*VM2;
			RM3=RM3+H*VM3;
			VM1=VM1+H*AM1;
			VM2=VM2+H*AM2;
			VM3=VM3+H*AM3;
			T=T+H;
		end
		RTM1=RT1-RM1;
		RTM2=RT2-RM2;
		RTM3=RT3-RM3;
		RTM=sqrt(RTM1^2+RTM2^2+RTM3^2);
		VTM1=VT1-VM1;
		VTM2=VT2-VM2;
		VTM3=VT3-VM3;
		VC=-(RTM1*VTM1+RTM2*VTM2+RTM3*VTM3)/RTM;
		TGO=RTM/VC;
		BETAD=XNT/VT;
		XNT1=XNT*sin(BETA);
		XNT2=XNT*cos(BETA);
		XNT3=0.;
		VT1=-VT*cos(BETA);
		VT2= VT*sin(BETA);
		VT3=0.;
		VM=sqrt(VM1*VM1+VM2*VM2+VM3*VM3);
		XNC1=((6.*RTM1+6.*VTM1*TGO)/TGO^2) +2.0*(VM1-VM1F)/(TGO)+XNT1;
		XNC2=((6.*RTM2+6.*VTM2*TGO)/TGO^2) +2.0*(VM2-VM2F)/(TGO)+XNT2;
		XNC3=((6.*RTM3+6.*VTM3*TGO)/TGO^2) +2.0*(VM3-VM3F)/(TGO)+XNT3;
		XNCG=sqrt(XNC1^2+XNC2^2+XNC3^2)/32.2;
		XNCDOTVM=(XNC1*VM1+XNC2*VM2+XNC3*VM3)/VM;
		AM1=XNC1-XNCDOTVM*VM1/VM;
		AM2=XNC2-XNCDOTVM*VM2/VM;
		AM3=XNC3-XNCDOTVM*VM3/VM;
		GAMYDEG=57.3*atan2(VM3,VM1);
		GAMPDEG=57.3*atan2(VM2,sqrt(VM1^2+VM3^2));
		FLAG=1;
	end
	FLAG=0;
	BETA=.5*(BETAOLD+BETA+H*BETAD);
	RT1=.5*(RT1OLD+RT1+H*VT1);
	RT2=.5*(RT2OLD+RT2+H*VT2);
	RT3=.5*(RT3OLD+RT3+H*VT3);
	RM1=.5*(RM1OLD+RM1+H*VM1);
	RM2=.5*(RM2OLD+RM2+H*VM2);
	RM3=.5*(RM3OLD+RM3+H*VM3);
	VM1=.5*(VM1OLD+VM1+H*AM1);
	VM2=.5*(VM2OLD+VM2+H*AM2);
	VM3=.5*(VM3OLD+VM3+H*AM3);
	S=S+H;
	if S>=.0999
		S=0.;
		n=n+1;
		RT1K=RT1/1000.;
		RT2K=RT2/1000.;
		RT3K=RT3/1000.;
		RM1K=RM1/1000.;
		RM2K=RM2/1000.;
		RM3K=RM3/1000.;
		ArrayT(n)=T;
		ArrayRM1K(n)=RM1K;
		ArrayRM2K(n)=RM2K;
		ArrayRM3K(n)=RM3K;
		ArrayRT1K(n)=RT1K;
		ArrayRT2K(n)=RT2K;
		ArrayRT3K(n)=RT3K;
		ArrayGAMPDEG(n)=GAMPDEG;
		ArrayGAMFPDEG(n)=GAMFPDEG;
		ArrayGAMYDEG(n)=GAMYDEG;
		ArrayGAMFYDEG(n)=GAMFYDEG;
		ArrayXNCG(n)=XNCG;
    end
end
n=n+1;
		RT1K=RT1/1000.;
		RT2K=RT2/1000.;
		RT3K=RT3/1000.;
		RM1K=RM1/1000.;
		RM2K=RM2/1000.;
		RM3K=RM3/1000.;
		ArrayT(n)=T;
		ArrayRM1K(n)=RM1K;
		ArrayRM2K(n)=RM2K;
		ArrayRM3K(n)=RM3K;
		ArrayRT1K(n)=RT1K;
		ArrayRT2K(n)=RT2K;
		ArrayRT3K(n)=RT3K;
		ArrayGAMPDEG(n)=GAMPDEG;
		ArrayGAMFPDEG(n)=GAMFPDEG;
		ArrayGAMYDEG(n)=GAMYDEG;
		ArrayGAMFYDEG(n)=GAMFYDEG;
		ArrayXNCG(n)=XNCG;
figure
plot(ArrayRM1K,ArrayRM2K,ArrayRT1K,ArrayRT2K),grid
xlabel('Downrange (kft)')
ylabel('Altitude (kft)')
figure
plot(ArrayRM1K,ArrayRM3K,ArrayRT1K,ArrayRT3K),grid
xlabel('Downrange (kft)')
ylabel('Crossrange (kft)')
figure
plot(ArrayT,ArrayGAMPDEG,ArrayT,ArrayGAMFPDEG),grid
xlabel('Time (Sec)')
ylabel('Pitch Reentry Angle (deg)')
figure
plot(ArrayT,ArrayGAMYDEG,ArrayT,ArrayGAMFYDEG),grid
xlabel('Time (Sec)')
ylabel('Yaw Reentry Angle (deg)')
figure
plot(ArrayT,ArrayXNCG),grid
xlabel('Time (Sec)')
ylabel('Commanded Acceleration (g)')
axis([00,10,00,60])
clc
output=[ArrayT',ArrayRM1K',ArrayRM2K',ArrayRM3K',ArrayRT1K',...
		ArrayRT2K',ArrayRT3K',ArrayGAMPDEG',ArrayGAMFPDEG',...
		ArrayGAMYDEG',ArrayGAMFYDEG'];
save datfil.txt output  -ascii
disp 'simulation finished'
RTM
		
 