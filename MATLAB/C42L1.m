% Runs very slowly because of small intergration interval
%Preallocation
clear
Z=zeros(size(1:1000));
I=zeros(size(1:50));
TF=zeros(size(1:50));
count=0;
TSW=1;
VC=5.*3280.;
XNTIC=161;
YIC=0.;
VM=3000.;
HEDEG=20.;
XNP=3.;
SIGNOISE=.0001;
TS=.01;
TAU=.2;
NOISE=1;
RUN=50;
% TYPE OF GUIDANCE (0=PN,1=OG,2=DG,3=HYBRID)
APN=1;
XLIM=322;
XNU=.5;
W=2;
% TYPE OF MANEUVER (1=POISSON,2=UNIF CONST,3=RANDOM SINE,4=RANDOM VS
ICONSTANT=2;
for TF=.1:.1:10.0,
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
        SUM=rand(1);
        PHASE=6.28*SUM;
		Y=YIC;
		YD=0;
		TS2=TS*TS;
		TS3=TS2*TS;
		TS4=TS3*TS;
		TS5=TS4*TS;
		PHIN=XNTIC*XNTIC/10;
		RTM=VC*TF;
		SIGPOS=RTM*SIGNOISE;
		SIGN2=SIGPOS^2;
		P11=SIGN2;
		P12=0.;
		P13=0.;
		P22=(VM*HEDEG/57.3)^2;
		P23=0.;
		P33=XNTIC*XNTIC;
		T=0.;
		H=.001;
		S=0.;
		YH=0.;
		YDH=0.;
		XNTH=0.;
		XNC=0.;
		XNL=0.;
		BETA=XNTIC;
		QFIRST=1;
		SIG=1./sqrt(2.*XNU);
		XNOISE=randn;
		if XNOISE>0
			XNTP=BETA;
		else
			XNTP=-BETA;
		end;
		DELT=9999.;
		TNOW=0.;
		while T <= (TF - 1e-5)
			if QFIRST==1
				XNOISE1=SIG*randn;
				XNOISE2=SIG*randn;
				DELT=XNOISE1^2+XNOISE2^2;
				QFIRST=0;
				TNOW=T;
			end;
			if T>= (DELT+TNOW)
				XNTP=-XNTP;
				QFIRST=1;
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
				end;
				if ICONSTANT==1
					XNTC=XNTP;
				elseif ICONSTANT==2
					if T<TSTART
						XNTC=0.;
					else
						XNTC=COEF*XNTIC;
					end;
				elseif ICONSTANT==3
					if T<TSTART
						XNTC=0.;
					else
						XNTC=XNTIC*sin(W*T+PHASE);
					end;
				else
					if T<TSTART
						XNTC=0.;
					else
						XNTC=COEF*XNTIC*sign(sin(W*(T-TSTART)));
					end;
				end;
				TGO=TF-T+.00001;
				RTM=VC*TGO;
				XLAM=Y/(VC*TGO);
				XLAMD=(RTM*YD+Y*VC)/(RTM^2);
				XNLD=(XNC-XNL)/TAU;
				YDD=XNTC-XNL;
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
      				M11=P11+TS*P12+.5*TS2*P13+TS*(P12+TS*P22+.5*TS2*P23);
      				M11=M11+.5*TS2*(P13+TS*P23+.5*TS2*P33)+TS5*PHIN/20.;
      				M12=P12+TS*P22+.5*TS2*P23+TS*(P13+TS*P23+.5*TS2*P33)...
					+TS4*PHIN/8.;
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
					XLAMNOISE=SIGNOISE*randn;
      				YSTAR=RTM*(XLAM+XLAMNOISE);
      				RES=YSTAR-YH-TS*YDH-.5*TS*TS*(XNTH-XNC);
      				YH=K1*RES+YH+TS*YDH+.5*TS*TS*(XNTH-XNC);      
      				YDH=K2*RES+YDH+TS*(XNTH-XNC);
      				XNTH=K3*RES+XNTH;
      				XLAMDH=(YH+YDH*TGO)/(VC*TGO*TGO);
					X=TGO/TAU;
					ZEM1H=YH+YDH*TGO-XNL*TAU*TAU*(exp(-X)+X-1.);
					ZEM2H=YH+YDH*TGO-XNL*TAU*TAU*(exp(-X)+X-1.)+...
					.5*XNTH*TGO*TGO;
      				if APN==0
         				XNC=XNP*(YH+YDH*TGO)/TGO^2;
      				elseif APN==1
         				X=TGO/TAU;
         				TOP=6.*X*X*(exp(-X)-1.+X);
         				BOT1=2*X*X*X+3.+6.*X-6.*X*X;
         				BOT2=-12.*X*exp(-X)-3.*exp(-2.*X);
         				XNPP=TOP/(.0001+BOT1+BOT2);
						XNC=XNPP*ZEM2H/TGO^2;
					elseif APN==2
						XNC=XLIM*sign(ZEM1H);
					else
						if TGO>TSW
							TOP=6.*X*X*(exp(-X)-1.+X);
							BOT1=2*X*X*X+3.+6.*X-6.*X*X;
							BOT2=-12.*X*exp(-X)-3.*exp(-2.*X);
							XNPP=TOP/(.0001+BOT1+BOT2);
							XNEW=XNPP*XNL*(exp(-X)+X-1.)/(X*X);
							XNC=XNPP*ZEM2H/TGO^2;
							if XNC>XLIM
								XNC=XLIM;
							end
							if XNC<-XLIM
								XNC=-XLIM;
							end
						else
							XNC=XLIM*sign(ZEM1H);
						end;
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
	Z2=0.;
 	for I=1:RUN,
       		Z1=(Z(I)-XMEAN)^2+Z1;
			Z2=Z(I)^2+Z2;
       		if I==1,
          		SIGMA=0;
				RMS=0.;
       		else
          		SIGMA=sqrt(Z1/(I-1));
				RMS=sqrt(Z2/(I-1));
       		end;
 	end;
	count=count+1;
	ArrayTF(count)=TF;
	ArraySIGMA(count)=SIGMA;
	ArrayXMEAN(count)=XMEAN;
	ArrayRMS(count)=RMS;
end;
figure
plot(ArrayTF',ArrayRMS'),grid
title('RMS miss for various flight times')
xlabel('Flight Time (S)')
ylabel('RMS MISS (Ft) ')
clc
output=[ArrayTF',ArrayRMS'];
save datfil.txt output -ascii
disp('Simulation Complete')
	
