clear
count=0;
VM=3000.;
VT=1000.;
XNT=0.;
RM1IC=0.;
RM2IC=1.;
RT1IC=40000.;
RT2IC=10000.;
HEDEG=0.;
XNP=10.;
TS=.1;
RT1=RT1IC;
RT2=RT2IC;
RM1=RM1IC;
RM2=RM2IC;
BETAT=0.;
VT1=-VT*cos(BETAT);
VT2=VT*sin(BETAT);
HE=HEDEG/57.3;
XNC=0.;
T=0.;
S=0.;
RTM1=RT1-RM1;
RTM2=RT2-RM2;
RTM=sqrt(RTM1^2+RTM2^2);
THETT=atan2(RT2,RT1);
VM1=VM*cos(THETT+HE);
VM2=VM*sin(THETT+HE);
VTM1=VT1-VM1;
VTM2=VT2-VM2;
VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
while VC >= 0        
  	if RTM < 1000
      		H=.0002;
   	else
      		H=.01;
   	end
   	BETATOLD=BETAT;
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
         		BETAT=BETAT+H*BETATD;
         		RT1=RT1+H*VT1;
         		RT2=RT2+H*VT2;
         		RM1=RM1+H*VM1;
         		RM2=RM2+H*VM2;
         		VM1=VM1+H*AM1;
         		VM2=VM2+H*AM2;
         		T=T+H;
      		end
      		THETT=atan2(RT2,RT1);
      		THETM=atan2(RM2,RM1);
      		RT=sqrt(RT1^2+RT2^2);
      		RM=sqrt(RM1^2+RM2^2);
      		RTM1=RT1-RM1;
      		RTM2=RT2-RM2;
      		RTM=sqrt(RTM1^2+RTM2^2);
      		VTM1=VT1-VM1;
      		VTM2=VT2-VM2;
      		VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
      		XLAM=atan2(RTM2,RTM1);
      		XLAMD=(RTM1*VTM2-RTM2*VTM1)/(RTM*RTM);
      		XNC=XNP*RM*(THETT-THETM);
      		AM1=-XNC*sin(XLAM);
      		AM2=XNC*cos(XLAM);
      		VT1=-VT*cos(BETAT);
      		VT2=VT*sin(BETAT);
      		BETATD=XNT/VT;
      		FLAG=1;
   	end
   	FLAG=0;
   	BETAT=.5*(BETATOLD+BETAT+H*BETATD);
   	RT1=.5*(RT1OLD+RT1+H*VT1);
   	RT2=.5*(RT2OLD+RT2+H*VT2);
   	RM1=.5*(RM1OLD+RM1+H*VM1);
   	RM2=.5*(RM2OLD+RM2+H*VM2);
   	VM1=.5*(VM1OLD+VM1+H*AM1);
   	VM2=.5*(VM2OLD+VM2+H*AM2);
   	S=S+H;
   	if S>=(TS - 1e-5)
   		S=0.;
      		RT1K=RT1/1000.;
      		RT2K=RT2/1000.;
      		RM1K=RM1/1000.;
      		RM2K=RM2/1000.;
      		count=count+1;
      		ArrayT(count)=T;
      		ArrayRT1K(count)=RT1K;
      		ArrayRT2K(count)=RT2K;
      		ArrayRM1K(count)=RM1K;
      		ArrayRM2K(count)=RM2K;
      		ArrayXNCG(count)=XNC/32.2;
   	end
end
RTM
figure
plot(ArrayRT1K,ArrayRT2K,ArrayRM1K,ArrayRM2K)
xlabel('Downrange (Kft)')
ylabel('Altitude (Kft) ')
figure
plot(ArrayT,ArrayXNCG),grid
xlabel('Time (S)')
ylabel('Acceleration (G) ')
clc
output=[ArrayT',ArrayRT1K',ArrayRT2K',ArrayRM1K',ArrayRM2K',ArrayXNCG'];
save datfil output -ascii
disp 'simulation finished'
