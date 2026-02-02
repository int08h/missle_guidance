n=0;
VM = 3000.;
VT = 1000.;
XNT = 0.;
HEDEG = -20.;
XNP = 4.;
RM1 = 0.;
RM2 = 10000.;
RT1 = 40000.;
RT2 = 10000.;
BETA=0.;
VT1=-VT*cos(BETA);
VT2=VT*sin(BETA);
HE=HEDEG/57.3;
T=0.;
S=0.;
RTM1=RT1-RM1;
RTM2=RT2-RM2;
RTM=sqrt(RTM1*RTM1+RTM2*RTM2);
XLAM=atan2(RTM2,RTM1); 
XLEAD=asin(VT*sin(BETA+XLAM)/VM);
THET=XLAM+XLEAD;
VM1=VM*cos(THET+HE);
VM2=VM*sin(THET+HE);
VTM1 = VT1 - VM1;
VTM2 = VT2 - VM2;
VC=-(RTM1*VTM1 + RTM2*VTM2)/RTM;
while VC >= 0 
  	if RTM < 1000
      		H=.0002;
   	else
      		H=.01;
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
      		RTM=sqrt(RTM1*RTM1+RTM2*RTM2);
      		VTM1=VT1-VM1;
      		VTM2=VT2-VM2;
      		VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
      		XLAM=atan2(RTM2,RTM1);
      		XLAMD=(RTM1*VTM2-RTM2*VTM1)/(RTM*RTM);
      		XNC=XNP*VC*XLAMD;
      		AM1=-XNC*sin(XLAM);
      		AM2=XNC*cos(XLAM);
      		VT1=-VT*cos(BETA);
      		VT2=VT*sin(BETA);
      		BETAD=XNT/VT;
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
		ArrayT(n)=T;
		ArrayRT1(n)=RT1;
		ArrayRT2(n)=RT2;
		ArrayRM1(n)=RM1;
		ArrayRM2(n)=RM2;
		ArrayXNCG(n)=XNC/32.2;
		ArrayRTM(n)=RTM;
	end
end
RTM
figure
plot(ArrayRT1,ArrayRT2,ArrayRM1,ArrayRM2),grid
title('Two-dimensional tactical missile-target engagement simulation')
xlabel('Downrange (Ft) ')
ylabel('Altitude (Ft)')
figure
plot(ArrayT,ArrayXNCG),grid
title('Two-dimensional tactical missile-target engagement simulation')
xlabel('Time (sec)')
ylabel('Acceleration of missle (G)')
clc
output=[ArrayT',ArrayRT1',ArrayRT2',ArrayRM1',ArrayRM2',ArrayXNCG',ArrayRTM' ];
save datfil output -ascii
disp '*** Simulation Complete'
