clear
XNP=3.;
TAU=1.;
XNT=193.2;
W=3.;
n=0;
for RT1IC=500:500:40000
   	VM=3000.;
      	VT=1000.;
      	RM1=0.;
      	RM2=0.;
      	RT1=RT1IC;
      	RT2=0.;      
      	BETA=0.;
      	VT1=-VT*cos(BETA);
      	VT2=VT*sin(BETA);
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
      	XLAMH=0.;
      	H=.01;
      	while VC>0.
         	if(RTM<1000.)
            		H=.0005;
         	end
            	BETAOLD=BETA;
            	RT1OLD=RT1;
            	RT2OLD=RT2;
            	RM1OLD=RM1;
            	RM2OLD=RM2;
            	VM1OLD=VM1;
            	VM2OLD=VM2;
            	XLAMHOLD=XLAMH;
            	STEP=1;
            	FLAG=0;
            	while STEP<=1
               	if FLAG==1
                  	STEP=2;
                  	BETA=BETA+H*BETAD;
                  	RT1=RT1+H*VT1;
                  	RT2=RT2+H*VT2;
                  	RM1=RM1+H*VM1;
                  	RM2=RM2+H*VM2;
                  	VM1=VM1+H*AM1;
                  	VM2=VM2+H*AM2;
                  	XLAMH=XLAMH+H*XLAMHD;
                  	T=T+H;
               	end
               	VT1=-VT*cos(BETA);
               	VT2=VT*sin(BETA);
               	BETAD=XNT*sin(W*T)/VT;
               	RTM1=RT1-RM1;
               	RTM2=RT2-RM2;
               	RTM=sqrt(RTM1^2+RTM2^2);
               	VTM1=VT1-VM1;
               	VTM2=VT2-VM2;
               	VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
               	XLAM=atan2(RTM2,RTM1);
               	XLAMHD=(XLAM-XLAMH)/TAU;
               	XNC=XNP*VC*XLAMHD;
               	AM1=-XNC*sin(XLAM);
               	AM2=XNC*cos(XLAM);
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
      	XLAMH=.5*(XLAMHOLD+XLAMH+H*XLAMHD);
   end
   if RTM2>0.
	RTMP=RTM;
   else
	RTMP=-RTM;
   end 
   n=n+1;
   ArrayT(n)=T;
   ArrayRTMP(n)=RTMP;
end
figure
plot(ArrayT,ArrayRTMP),grid
xlabel('Flight Time (Sec)')
ylabel('Miss (Ft)')
clc
output=[ArrayT',ArrayRTMP'];
save datfil.txt output  -ascii
disp 'simulation finished'

