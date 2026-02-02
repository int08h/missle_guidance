clear
%Preallocation
Z=zeros(size(1:1000));
I=zeros(size(1:50));
TF=zeros(size(1:50));
count=0;
VC=4000.;
XNT=96.6;
YIC=0.;
VM=3000.;
HEDEG=20.;
XNP=3.;
SIGNOISE=.001;
TS=.1;
TAU=.5;
NOISE=1;
RUN=50;
APN=0;
XLIM=999999;
for TF=.5:.5:10.0,
	Z1=0;
	for I=1:RUN
		Y=YIC;
		YD=0;
		TS2=TS*TS;
		TS3=TS2*TS;
		TS4=TS3*TS;
		TS5=TS4*TS;
		PHIN=XNT*XNT/TF;
		RTM=VC*TF;
		SIGPOS=RTM*SIGNOISE;
		SIGN2=SIGPOS^2;
		P11=SIGN2;
		P12=0.;
		P13=0.;
		P22=(VM*HEDEG/57.3)^2;
		P23=0.;
		P33=XNT*XNT;
		T=0.;
		H=.01;
		S=0.;
		YH=0.;
		YDH=0.;
		XNTH=0.;
		XNC=0.;
		XNL=0.;
		while T <= (TF - 1e-5)
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
      				end;
      				TGO=TF-T+.00001;
      				RTM=VC*TGO;
      				XLAM=Y/(VC*TGO);
      				XLAMD=(RTM*YD+Y*VC)/(RTM^2);
      				XNLD=(XNC-XNL)/TAU;
      				YDD=XNT-XNL;
      				FLAG=1;
   			end;
   			FLAG=0;
 			Y=.5*(YOLD+Y+H*YD);
 			YD=.5*(YDOLD+YD+H*YDD);
 			XNL=.5*(XNLOLD+XNL+H*XNLD);
			S=S+H;
   			if S>=(TS - 1e-5)
      				S=0.;
      				TGO=TF-T+.000001;
      				RTM=VC*TGO;
      				SIGPOS=RTM*SIGNOISE;
      				SIGN2=SIGPOS^2;
      				M11=P11+TS*P12+.5*TS2*P13+TS*(P12+TS*P22+...
							.5*TS2*P23);
      				M11=M11+.5*TS2*(P13+TS*P23+.5*TS2*P33)...
						+TS5*PHIN/20.;
      				M12=P12+TS*P22+.5*TS2*P23+TS*(P13+TS*P23+...
						.5*TS2*P33)+TS4*PHIN/8.;
      				M13=P13+TS*P23+.5*TS2*P33+PHIN*TS3/6.;
      				M22=P22+TS*P23+TS*(P23+TS*P33)+PHIN*TS3/3.;
      				M23=P23+TS*P33+.5*TS2*PHIN;
      				M33=P33+PHIN*TS;
      				K1=M11/(M11+SIGN2);
      				K2=M12/(M11+SIGN2);
      				K3=M13/(M11+SIGN2);
      				P11=(1.-K1)*M11;
      				P12=(1.-K1)*M12;
      				P13=(1.-K1)*M13;
      				P22=-K2*M12+M22;
      				P23=-K2*M13+M23;
      				P33=-K3*M13+M33;
					if NOISE==1,
							XLAMNOISE=SIGNOISE*randn;
					else
                			XLAMNOISE=0;
					end;
      				YSTAR=RTM*(XLAM+XLAMNOISE);
      				RES=YSTAR-YH-TS*YDH-.5*TS*TS*(XNTH-XNC);
      				YH=K1*RES+YH+TS*YDH+.5*TS*TS*(XNTH-XNC);      
      				YDH=K2*RES+YDH+TS*(XNTH-XNC);
      				XNTH=K3*RES+XNTH;
      				XLAMDH=(YH+YDH*TGO)/(VC*TGO*TGO);
      				if APN==0
         				XNC=XNP*VC*XLAMDH;
      				elseif APN==1
         				XNC=XNP*VC*XLAMDH+APN*.5*XNP*XNTH;
      				else
         				X=TGO/TAU;
         				TOP=6.*X*X*(exp(-X)-1.+X);
         				BOT1=2*X*X*X+3.+6.*X-6.*X*X;
         				BOT2=-12.*X*exp(-X)-3.*exp(-2.*X);
         				XNPP=TOP/(.0001+BOT1+BOT2);
         				XNEW=XNPP*XNL*(EXP(-X)+X-1.)/(X*X);
					XNC=XNPP*VC*XLAMDH+.5*XNPP*XNTH-XNEW;
      				end;
      				if XNC>XLIM
      					XNC=XLIM;
      				elseif XNC<-XLIM
      					XNC=-XLIM;
      				end;
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
ylabel('Mean Miss (Ft) ')
axis([00,10,00,60])
clc
output=[ArrayTF',ArraySIGMA',ArrayXMEAN'];
save datfil.txt output -ascii
disp('Simulation Complete')
	
