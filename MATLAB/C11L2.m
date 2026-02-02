clear
count=0;
APN=0;
XNP=3.;
RT1=0.;
RT2=200000.;
RM1=170000.;
RM2=0.;
VT=6000.;
RT2DES=50000.;
GAMTDEG=45.;
BETA=500.;
BETEST=500.;
XNCLIMG=7.;
XNCLIM=XNCLIMG*32.2;
VT1=VT*cos(GAMTDEG/57.3);
VT2=-VT*sin(GAMTDEG/57.3);
[RT1F,RT2F,TFDES]=initialpz(RT2DES,RT1,RT2,VT1,VT2,BETEST);
RTM1F=RT1F-RM1;
RTM2F=RT2F-RM2;
GAMMDEG=57.3*atan2(RTM2F,RTM1F);
RTMF=sqrt(RTM1F^2+RTM2F^2);
VM=RTMF/TFDES;
VM1=VM*cos(GAMMDEG/57.3);
VM2=VM*sin(GAMMDEG/57.3);
RTM1=RT1-RM1;
RTM2=RT2-RM2;
RTM=sqrt(RTM1^2+RTM2^2);
VTM1=VT1-VM1;
VTM2=VT2-VM2;
VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
T=0.;
H=.01;
S=0.;
XNC=0.;
ZEMPLOS=0.;
ZEM1=0.;
ZEM2=0.;
while VC>=0.
 	if RTM<1000.
		H=.0002;
	else
		H=.01;
	end
	RT1OLD=RT1;
	RT2OLD=RT2;
	VT1OLD=VT1;
	VT2OLD=VT2;
	RM1OLD=RM1;
	RM2OLD=RM2;
	VM1OLD=VM1;
	VM2OLD=VM2;
 	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
			RT1=RT1+H*VT1;
			RT2=RT2+H*VT2;
			VT1=VT1+H*AT1;
			VT2=VT2+H*AT2;
			RM1=RM1+H*VM1;
			RM2=RM2+H*VM2;
			VM1=VM1+H*AM1;
			VM2=VM2+H*AM2;
 			T=T+H;
 			STEP=2;
 		end
 		if RT2<=30000.
			RHO=.002378*exp(-RT2/30000.);
		else
			RHO=.0034*exp(-RT2/22000.);
		end
		VT=sqrt(VT1^2+VT2^2);
		Q=.5*RHO*VT^2;
		GAMT=atan2(-VT2,VT1);
		AT1=-32.2*Q*cos(GAMT)/BETA;
		AT2=-32.2+32.2*Q*sin(GAMT)/BETA;
		RTM1=RT1-RM1;
		RTM2=RT2-RM2;
		RTM=sqrt(RTM1^2+RTM2^2);
		VTM1=VT1-VM1;
		VTM2=VT2-VM2;
		VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
		XLAM=atan2(RTM2,RTM1);
		XLAMD=(RTM1*VTM2-RTM2*VTM1)/(RTM*RTM);
		ATPLOS=-AT1*sin(XLAM)+AT2*cos(XLAM);
		if APN==1
			XNC=XNP*VC*XLAMD+.5*XNP*ATPLOS;
		else
			XNC=XNP*VC*XLAMD;
		end
		if XNC>XNCLIM
			XNC=XNCLIM;
		end
		if XNC<-XNCLIM
			XNC=-XNCLIM;
		end
		AM1=-XNC*sin(XLAM);
		AM2=XNC*cos(XLAM);
		FLAG=1;
	end;
	FLAG=0;
	RT1=.5*(RT1OLD+RT1+H*VT1);
	RT2=.5*(RT2OLD+RT2+H*VT2);
	VT1=.5*(VT1OLD+VT1+H*AT1);
	VT2=.5*(VT2OLD+VT2+H*AT2);
	RM1=.5*(RM1OLD+RM1+H*VM1);
	RM2=.5*(RM2OLD+RM2+H*VM2);
	VM1=.5*(VM1OLD+VM1+H*AM1);
	VM2=.5*(VM2OLD+VM2+H*AM2);
	S=S+H;
 	if S>=.09999
		S=0.;
		ATG=sqrt(AT1^2+AT2^2)/32.2;
		RT1K=RT1/1000.;
		RT2K=RT2/1000.;
		RM1K=RM1/1000.;
		RM2K=RM2/1000.;
		XNCG=XNC/32.2;
		ATPLOSG=ATPLOS/32.2;
		count=count+1;
		ArrayT(count)=T;
		ArrayRT1K(count)=RT1K;
		ArrayRT2K(count)=RT2K;
		ArrayRM1K(count)=RM1K;
		ArrayRM2K(count)=RM2K;
		ArrayATG(count)=ATG;
		ArrayATPLOSG(count)=ATPLOSG;
		ArrayXNCG(count)=XNCG;
 	end
 end
 RTM
figure
plot(ArrayRT1K,ArrayRT2K,ArrayRM1K,ArrayRM2K),grid
xlabel('Downrange (Kft)')
ylabel('Altitude (Kft) ')
figure
plot(ArrayT,ArrayATG,ArrayT,ArrayATPLOSG,ArrayT,ArrayXNCG),grid
xlabel('Time (Sec)')
ylabel('Acceleration (G) ')
clc
output=[ArrayT',ArrayRT1K',ArrayRT2K',ArrayRM1K',ArrayRM2K',ArrayATG',ArrayATPLOSG',ArrayXNCG'];
save datfil.txt output -ascii
disp 'simulation finished'
