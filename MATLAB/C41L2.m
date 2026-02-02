%Preallocation
clear
Z=zeros(size(1:1000));
I=zeros(size(1:50));
TF=zeros(size(1:50));
count=0;		
VC=5.*3280.;
XNTIC=161.;
YIC=0.;
VM=3000.;
HEDEG=20.;
XNP=3.;
SIGNOISE=10.;
TS=.01;
TAU=.2;
RUN=100;
AMAXG=99999999.;
PZ1=.0001;
AMAX=AMAXG*32.2;
PHIN=SIGNOISE*SIGNOISE*TS;
for TF=.2:.2:10.0,
	Z1=0.;
	for JJ=1:RUN
		SUM=rand(1);
		TSTART=TF*SUM;
		PZ=rand(1);
		PZ=PZ-.5;
		if PZ > 0
			COEF=1;
		else
			COEF=-1;
		end;
		Y=0.;
		YD=0.;
		TS2=TS*TS;
		TS3=TS2*TS;
		TS4=TS3*TS;
		TS5=TS4*TS;
		PHIS=XNTIC*XNTIC/TF;
		RTM=VC*TF;
		SIGPOS=SIGNOISE;
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
				if T<TSTART
					XNTC=0.;
				else
					XNTC=COEF*XNTIC;
				end;
 				TGO=TF-T+.00001;
				RTM=VC*TGO;
				XLAM=Y/(VC*TGO);
				XNLD=(XNC-XNL)/TAU;
				XNT=XNTC;
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
				SIGPOS=SIGNOISE;
				SIGN2=SIGPOS^2;
				M11=P11+TS*P12+.5*TS2*P13+TS*(P12+TS*P22+.5*TS2*P23);
				M11=M11+.5*TS2*(P13+TS*P23+.5*TS2*P33)+TS5*PHIS/20.;
				M12=P12+TS*P22+.5*TS2*P23+TS*(P13+TS*P23+.5*TS2*P33);
     				M12=M12+TS4*PHIS/8.;
				M13=P13+TS*P23+.5*TS2*P33+PHIS*TS3/6.;
				M22=P22+TS*P23+TS*(P23+TS*P33)+PHIS*TS3/3.;
				M23=P23+TS*P33+.5*TS2*PHIS;
				M33=P33+PHIS*TS;
				K1=M11/(M11+SIGN2);
				K2=M12/(M11+SIGN2);
				K3=M13/(M11+SIGN2);
				P11=(1.-K1)*M11;
				P12=(1.-K1)*M12;
				P13=(1.-K1)*M13;
				P22=-K2*M12+M22;
				P23=-K2*M13+M23;
				P33=-K3*M13+M33;
				YNOISE=SIGNOISE*randn;
				YSTAR=Y+YNOISE;
				RES=YSTAR-YH-TS*YDH-.5*TS*TS*(XNTH-XNL);
				YH=K1*RES+YH+TS*YDH+.5*TS*TS*(XNTH-XNL);
				YDH=K2*RES+YDH+TS*(XNTH-XNL);
				XNTH=K3*RES+XNTH;
				X=TGO/TAU;
				ZEM2H=YH+YDH*TGO-XNL*TAU*TAU*(exp(-X)+X-1.)...
                    +.5*XNTH*TGO*TGO;
				TOP=6.*X*X*(exp(-X)-1.+X);
				BOT1=2*X*X*X+3.+6.*X-6.*X*X;
				BOT2=-12.*X*exp(-X)-3.*exp(-2.*X);
				XNPP=TOP/(PZ1+BOT1+BOT2);
				XNEW=XNPP*XNL*(exp(-X)+X-1.)/(X*X);
				XNC=XNPP*ZEM2H/TGO^2;
				if XNC>AMAX
					XNC=AMAX;
				end;
				if XNC<-AMAX
					XNC=-AMAX;
				end;
			end;

 		end;
		SP11=sqrt(P11);
		Z(JJ)=Y;
		Z1=Z(JJ)+Z1;
		XMEAN=Z1/JJ;
	end;
 	SIGMA=0.;
	Z1=0.;
	Z2=0.;
	for JJ=1:RUN,
		Z1=(Z(JJ)-XMEAN)^2+Z1;
		Z2=Z(JJ)^2+Z2;
		if JJ==1,
			SIGMA=0.;
			RMS=0.;
		else
			SIGMA=sqrt(Z1/(JJ-1));
			RMS=sqrt(Z2/(JJ-1));
		end;
	end;
 	FORM=sqrt(2.*(PHIS^.16667)*(PHIN^.8333));
 	count=count+1;
	ArrayTF(count)=TF;
	ArrayRMS(count)=RMS;
	ArraySP11(count)=SP11;
	ArrayFORM(count)=FORM;
end;
figure
plot(ArrayTF',ArrayRMS',ArrayTF',ArrayFORM',ArrayTF',ArraySP11'),grid
title('RMS miss for various flight times')
xlabel('Flight Time (S)')
ylabel('RMS MISS (Ft) ')
clc
output=[ArrayTF',ArrayRMS',ArraySP11',ArrayFORM'];
save datfil.txt output -ascii
disp('Simulation Complete')
