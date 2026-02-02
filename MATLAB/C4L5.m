% Preallocation
clear
Z=zeros(1,1000);
TF=zeros(1,10);
I=zeros(1,100);
ArrayTF=zeros(1,10);
ArraySIGMA=zeros(1,10);
ArrayXMEAN=zeros(1,10);
count=0;
VC=4000;
XNT=96.6;
VM=3000;
XNP=3;
TAU=1;
RUN=50
for TF=1:10
	Z1=0;
	for I=1:RUN
		SUM=rand(1);
		TSTART=TF*SUM;
		PZ=rand(1);
		PZ=PZ-.5;
		if PZ > 0
			COEF=1;
		else
			COEF=-1;
		end;
      	Y=0;
      	YD=0;
      	T=0;
      	H=.01;
      	S=0;
      	XNC=0;
      	XNL=0;
      	while T<=(TF - 1e-5)
         	if T < TSTART 
            	XNT=0;
         	else
            	XNT=COEF*96.6;
         	end;	
        	YOLD=Y;
        	YDOLD=YD;
        	XNLOLD=XNL;
        	STEP=1;
        	FLAG=0;
        	while STEP <=1
           	if FLAG==1
              	Y=Y+H*YD;
              	YD=YD+H*YDD;
              	XNL=XNL+H*XNLD;
              	T=T+H;
              	STEP=2;
           end
           TGO=TF-T+.00001;
           RTM=VC*TGO;
           XLAMD=(RTM*YD+Y*VC)/(RTM^2);
           XNC=XNP*VC*XLAMD;
           XNLD=(XNC-XNL)/TAU;
           YDD=XNT-XNL;
           FLAG=1;
        end;
        FLAG=0;
        Y=.5*(YOLD+Y+H*YD);
        YD=.5*(YDOLD+YD+H*YDD);
        XNL=.5*(XNLOLD+XNL+H*XNLD);
        S=S+H;
     end;
     Z(I)=Y;
     Z1=Z(I)+Z1;
     XMEAN=Z1/I;
	end;
   	SIGMA=0;
   	Z1=0;
   	for I=1:RUN
      	Z1=(Z(I)-XMEAN)^2+Z1;
      	if I==1 
         	SIGMA=0;
      	else
         	SIGMA=sqrt(Z1/(I-1));
      	end
   	end;
   	count=count+1;
   	ArrayTF(count)=TF;
   	ArraySIGMA(count)=SIGMA;
   	ArrayXMEAN(count)=XMEAN;
end; 
figure
plot(ArrayTF,ArraySIGMA,'r+') 
title('Shaping filter Monte Carlo results')
xlabel('Time')
ylabel('Standard Deviation / Mean')
axis([00,10,00,30])
clc
output=[ArrayTF',ArraySIGMA',ArrayXMEAN'];
save datfil.txt output  -ascii
disp 'simulation finished'
