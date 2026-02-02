clear
%Preallocation
Z=zeros(size(1:1000));
I=zeros(size(1:50));
TF=zeros(size(1:50));
count=0;
VC=4000;
XNT=96.6;
YIC=0;
VM=3000;
HEDEG=0;
BETA=.8;
XNP=3;
SIGNOISE=.001;
TS=.1;
RUN=50;
NOISE=1;
for TF=.5:.5:10.0,
	Z1=0;
	for I=1:RUN
		Y=YIC;
		YD=-VM*HEDEG/57.3;
		YDIC=YD;
		T=0.;
		H=.01;
		S=0.;
		GFILTER=1.-BETA^2;
		HFILTER=(1.-BETA)^2;
		XLAMH=0.;
		XLAMDH=0.;
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
				end;
    				TGO=TF-T+.00001;
				RTM=VC*TGO;
				XLAM=Y/(VC*TGO);
				XLAMD=(RTM*YD+Y*VC)/(RTM^2);
				YDD=XNT-XNC; 	
				FLAG=1;
			end;
			FLAG=0;
			Y=.5*(YOLD+Y+H*YD);
 			YD=.5*(YDOLD+YD+H*YDD);
			S=S+H;
 			if S>=(TS - 1e-5)
             			S=0.;
             			if NOISE==1,
                			XLAMNOISE=gaussc7(SIGNOISE);
             			else
                			XLAMNOISE=0;
             			end;
             			RES=XLAM-(XLAMH+TS*XLAMDH)+XLAMNOISE;
             			XLAMH=GFILTER*RES+XLAMH+TS*XLAMDH;
             			XLAMDH=HFILTER*RES/TS+XLAMDH;
             			XNC=XNP*VC*XLAMDH;
 			end;
		end;
		Z(I)=Y;
		Z1=Z(I)+Z1;
		XMEAN=Z1/I;
	end;
	SIGMA=0;
	Z1=0;
 	for I=1:RUN,
       		Z1=(Z(I)-XMEAN)^2+Z1;
       		if I==1,
          		SIGMA=0;
       		else
          		SIGMA=sqrt(Z1/(I-1));
       		end;
 	end;
	count=count+1;
	ArrayTF(count)=TF;
	ArraySIGMA(count)=SIGMA;
	ArrayXMEAN(count)=XMEAN;
end;
figure
plot(ArrayTF',ArraySIGMA'),grid
title('Standard deviation of miss for various flight times')
xlabel('Flight Time (S)')
ylabel('Noise Miss Standard Deviation (Ft) ')
axis([00,10,00,4])
figure
plot(ArrayTF',ArrayXMEAN'),grid
title('Mean of miss for various flight times')
xlabel('Flight Time (S)')
ylabel('Noise Miss Standard Deviation (Ft) ')
axis([00,10,00,60])
clc
output=[ArrayTF',ArraySIGMA',ArrayXMEAN'];
save datfil.txt output -ascii
disp('Simulation Complete')


