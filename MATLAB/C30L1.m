clear
n=0;
TAU=.5;
APN=0;
VC=9000.;
XNT=96.6;
XNTREAL=96.6;
XNTMAX=96.6;
W=2.;
YIC=0.;
VM=3000.;
HEDEG=0.;
HEDEGFIL=20.;
XNP=3.;
SIGRIN=.001;
TS=.01;
TF=10.;
Y=YIC;
YD=-VM*HEDEG/57.3;
YDIC=YD;
TS2=TS*TS;
TS3=TS2*TS;
TS4=TS3*TS;
TS5=TS4*TS;
PHIN=XNTMAX*XNTMAX/TF;
RTM=VC*TF;
SIGNOISE=SIGRIN;
SIGPOS=RTM*SIGNOISE;
SIGN2=SIGPOS^2;
P11=SIGN2;
P12=0.;
P13=0.;
P22=(VM*HEDEGFIL/57.3)^2;
P23=0.;
P33=XNTMAX*XNTMAX;
T=0.;
H=.001;
S=0.;
YH=0.;
YDH=0.;
XNTH=0.;
XNC=0.;
XNL=0.;
while T<=TF
 	YOLD=Y;
	YDOLD=YD;
	XNLOLD=XNL;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
         		STEP=2;
 			Y=Y+H*YD;
 			YD=YD+H*YDD;
			XNL=XNL+H*XNLD;
			T=T+H;
		end
		XNT=XNTREAL*sin(W*T);
 		TGO=TF-T+.00001;
		RTM=VC*TGO;
		XLAM=Y/(VC*TGO);
		XLAMD=(RTM*YD+Y*VC)/(RTM^2);
		XNLD=(XNC-XNL)/TAU;
		YDD=XNT-XNL;
		FLAG=1;
	end
	FLAG=0;
 	Y=.5*(YOLD+Y+H*YD);
 	YD=.5*(YDOLD+YD+H*YDD);
	XNL=.5*(XNLOLD+XNL+H*XNLD);
	S=S+H;
	if S>=(TS-.0001)
		S=0.;
		TGO=TF-T+.000001;
		RTM=VC*TGO;
		SIGNOISE=SIGRIN;
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
		RES=YSTAR-YH-TS*YDH-.5*TS*TS*(XNTH-XNL);
		YH=K1*RES+YH+TS*YDH+.5*TS*TS*(XNTH-XNL);
		YDH=K2*RES+YDH+TS*(XNTH-XNL);
		XNTH=K3*RES+XNTH;
		XLAMDH=(YH+YDH*TGO)/(VC*TGO*TGO);
		if APN==0
			XNC=XNP*(YH+YDH*TGO)/(TGO*TGO);
		elseif APN==1
			XNC=XNP*(YH+YDH*TGO+.5*XNTH*TGO*TGO)/(TGO*TGO);
		else
			XS=TGO/TAU;
			TOP=6.*XS*XS*(exp(-XS)-1.+XS);
			BOT1=2*XS*XS*XS+3.+6.*XS-6.*XS*XS;
			BOT2=-12.*XS*exp(-XS)-3.*exp(-2.*XS);
			XNPP=TOP/(.0001+BOT1+BOT2);
			C1=XNPP/(TGO*TGO);
			C2=XNPP/TGO;
			C3=.5*XNPP;
			C4=-XNPP*(exp(-XS)+XS-1.)/(XS*XS);
			XNC=C1*YH+C2*YDH+C3*XNTH+C4*XNL;
		end
		n=n+1;
		XNTG=XNT/32.2;
		XNTHG=XNTH/32.2;
		ArrayT(n)=T;
		ArrayXNTG(n)=XNTG;
		ArrayXNTHG(n)=XNTHG;
		ArrayY(n)=Y;
		ArrayYSTAR(n)=YSTAR;
	end
end
figure
plot(ArrayT,ArrayXNTG,ArrayT,ArrayXNTHG),grid
title('Acceleration Estimate')
xlabel('Time (Sec) ')
ylabel('Acceleration (G)')
figure
plot(ArrayT,ArrayY,ArrayT,ArrayYSTAR),grid
title('Measurement and Signal')
xlabel('Time (Sec) ')
ylabel('Y (Ft)')
clc
output=[ArrayT',ArrayXNTG',ArrayXNTHG',ArrayY',ArrayYSTAR'];
save datfil.txt output /ascii
disp '*** Simulation Complete'
