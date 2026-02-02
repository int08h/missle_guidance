clear
n=0;
XNP=3.;
DISPLACE=200.;
THOM=1.;
TAU=1.;
VM=3000.;
VT=1000.;
for THOM = 0.1:0.1:5
	RM1=0.;
	RM2=1000.;
	RT1=20000.;
	RT2=1000.;
	QSWITCH=0;	
	VT1=-VT;
	VT2=0.;
	T=0.;
	S=0.;
	RTM1=RT1-RM1;
	RTM2=RT2-RM2;
	RTM=sqrt(RTM1^2+RTM2^2);
	XLAM=atan2(RTM2,RTM1);
	VM1=VM;
	VM2=0.;
	VTM1=VT1-VM1;
	VTM2=VT2-VM2;
	VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
	TGO=RTM/VC;
	H=.001;
	D=0.;
	while VC>0.
		TGO=RTM/VC;
		if(TGO<.3)
				H=.00001;
		end
		if(TGO<=THOM & QSWITCH==0)
			QSWITCH=1;
			RT2=RT2+DISPLACE;
			RTM2=RT2-RM2;
			XLAM=atan2(RTM2,RTM1);
			D=XLAM;
		end	
		RT1OLD=RT1;
		RT2OLD=RT2;
		RM1OLD=RM1;
		RM2OLD=RM2;
		VM1OLD=VM1;
		VM2OLD=VM2;
		DOLD=D;
		STEP=1;
		FLAG=0;
		while STEP<=1
      		if FLAG==1
         		STEP=2;
         		RT1=RT1+H*VT1;
         		RT2=RT2+H*VT2;
         		RM1=RM1+H*VM1;
         		RM2=RM2+H*VM2;
         		VM1=VM1+H*AM1;
         		VM2=VM2+H*AM2;
				D=D+H*DD;
         		T=T+H;
      		end
      		RTM1=RT1-RM1;
      		RTM2=RT2-RM2;
      		RTM=sqrt(RTM1^2+RTM2^2);
      		VTM1=VT1-VM1;
      		VTM2=VT2-VM2;
      		VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
      		XLAM=atan2(RTM2,RTM1);
			EPS=XLAM-D;
			DD=(XLAM-D)/TAU;
			XNC=XNP*VC*DD;
      		AM1=-XNC*sin(XLAM);
      		AM2=XNC*cos(XLAM);
      		FLAG=1;
		end
		FLAG=0;
		RT1=.5*(RT1OLD+RT1+H*VT1);
		RT2=.5*(RT2OLD+RT2+H*VT2);
		RM1=.5*(RM1OLD+RM1+H*VM1);
		RM2=.5*(RM2OLD+RM2+H*VM2);
		VM1=.5*(VM1OLD+VM1+H*AM1);
		VM2=.5*(VM2OLD+VM2+H*AM2);
		D=.5*(DOLD+D+H*DD);
	end
	if RTM2>0
		RTMP=RTM;
	else
		RTMP=-RTM;
	end
	n=n+1;
	ArrayTHOM(n)=THOM;
	ArrayRTMP(n)=RTMP;
end
figure
plot(ArrayTHOM,ArrayRTMP),grid
xlabel('Homing Time (s)')
ylabel('Miss (Ft)')
clc
output=[ArrayTHOM',ArrayRTMP'];
save datfil.txt output  -ascii
disp 'simulation finished'


