clear
count=0;
TAU=.5;
ORDER=4;
VC=9000.;
XLIM=322.;
W1=1.;
W2=2.;
W3=4.;
WREAL=2.;
XNT=96.6;
XNTREAL=96.6;
TS=.01;
YIC=0.;
VM=3000.;
HEDEG=0.;
HEDEGFIL=20.;
XNP=3.;
SIGRIN=.0001;
TF=10.;;
PHASE=0./57.3;
X1=W1*TS;
X2=W2*TS;
X3=W3*TS;
Y=YIC;
YD=-VM*HEDEG/57.3;
PHIS1=W1*W1*XNT*XNT/TF;
PHIS2=W2*W2*XNT*XNT/TF;
PHIS3=W3*W3*XNT*XNT/TF;
RTM=VC*TF;
SIGNOISE=SIGRIN;
SIGPOS=RTM*SIGNOISE;
SIGN2=SIGPOS^2;
PHI1=zeros([4,4]);
P1=zeros([4,4]);
Q1=zeros([4,4]);
IDNP=eye(4);
PHI2=zeros([4,4]);
P2=zeros([4,4]);
Q2=zeros([4,4]);
PHI3=zeros([4,4]);
P3=zeros([4,4]);
Q3=zeros([4,4]);

PHI1(1,1)=1;
PHI1(1,2)=TS;
PHI1(1,3)=(1-cos(X1))/(W1*W1);
PHI1(1,4)=(X1-sin(X1))/(W1*W1*W1);
PHI1(2,2)=1;
PHI1(2,3)=sin(X1)/W1;
PHI1(2,4)=(1-cos(X1))/(W1*W1);
PHI1(3,3)=cos(X1);
PHI1(3,4)=sin(X1)/W1;
PHI1(4,3)=-W1*sin(X1);
PHI1(4,4)=cos(X1);
	
PHI2(1,1)=1;
PHI2(1,2)=TS;
PHI2(1,3)=(1-cos(X2))/(W2*W2);
PHI2(1,4)=(X2-sin(X2))/(W2*W2*W2);
PHI2(2,2)=1;
PHI2(2,3)=sin(X2)/W2;
PHI2(2,4)=(1-cos(X2))/(W2*W2);
PHI2(3,3)=cos(X2);
PHI2(3,4)=sin(X2)/W2;
PHI2(4,3)=-W2*sin(X2);
PHI2(4,4)=cos(X2);
	
PHI3(1,1)=1;
PHI3(1,2)=TS;
PHI3(1,3)=(1-cos(X3))/(W3*W3);
PHI3(1,4)=(X3-sin(X3))/(W3*W3*W3);
PHI3(2,2)=1;
PHI3(2,3)=sin(X3)/W3;
PHI3(2,4)=(1-cos(X3))/(W3*W3);
PHI3(3,3)=cos(X3);
PHI3(3,4)=sin(X3)/W3;
PHI3(4,3)=-W3*sin(X3);
PHI3(4,4)=cos(X3);
	
Q1(1,1)=PHIS1*(.333*X1^3-2*sin(X1)+2*X1*cos(X1)+.5*X1-...
			   .25*sin(2*X1))/(W1^5);
Q1(1,2)=PHIS1*(.5*X1*X1-X1*sin(X1)+.5*sin(X1)*...
			   sin(X1))/(W1^4);
Q1(2,1)=Q1(1,2);
Q1(1,3)=PHIS1*(sin(X1)-X1*cos(X1)-.5*X1+...
			   .25*sin(2*X1))/(W1^3);
Q1(3,1)=Q1(1,3);
Q1(1,4)=PHIS1*(cos(X1)+X1*sin(X1)-.5*sin(X1)*...
			   sin(X1)-1)/(W1*W1);
Q1(4,1)=Q1(1,4);
Q1(2,2)=PHIS1*(1.5*X1-2*sin(X1)+.25*sin(2*X1))/(W1^3);
Q1(2,3)=PHIS1*(1-cos(X1)-.5*sin(X1)*sin(X1))/(W1*W1);
Q1(3,2)=Q1(2,3);
Q1(2,4)=PHIS1*(sin(X1)-.5*X1-.25*sin(2*X1))/W1;
Q1(4,2)=Q1(2,4);
Q1(3,3)=PHIS1*(.5*X1-.25*sin(2*X1))/W1;
Q1(3,4)=.5*PHIS1*sin(X1)*sin(X1);
Q1(4,3)=Q1(3,4);
Q1(4,4)=W1*PHIS1*(.5*X1+.25*sin(2*X1));
	
Q2(1,1)=PHIS2*(.333*X2^3-2*sin(X2)+2*X2*cos(X2)...
			   +.5*X2-.25*sin(2*X2))/(W2^5);
Q2(1,2)=PHIS2*(.5*X2*X2-X2*sin(X2)+.5*sin(X2)...
			   *sin(X2))/(W2^4);
Q2(2,1)=Q2(1,2);
Q2(1,3)=PHIS2*(sin(X2)-X2*cos(X2)-.5*X2+.25*...
			   sin(2*X2))/(W2^3);
Q2(3,1)=Q2(1,3);
Q2(1,4)=PHIS2*(cos(X2)+X2*sin(X2)-.5*sin(X2)*...
			   sin(X2)-1)/(W2*W2);
Q2(4,1)=Q2(1,4);
Q2(2,2)=PHIS2*(1.5*X2-2*sin(X2)+.25*sin(2*X2))/(W2^3);
Q2(2,3)=PHIS2*(1-cos(X2)-.5*sin(X2)*sin(X2))/(W2*W2);
Q2(3,2)=Q2(2,3);
Q2(2,4)=PHIS2*(sin(X2)-.5*X2-.25*sin(2*X2))/W2;
Q2(4,2)=Q2(2,4);
Q2(3,3)=PHIS2*(.5*X2-.25*sin(2*X2))/W2;
Q2(3,4)=.5*PHIS2*sin(X2)*sin(X2);
Q2(4,3)=Q2(3,4);
Q2(4,4)=W2*PHIS2*(.5*X2+.25*sin(2*X2));
	
Q3(1,1)=PHIS3*(.333*X3^3-2*sin(X3)+2*X3*cos(X3)+...
			   .5*X3-.25*sin(2*X3))/(W3^5);
Q3(1,2)=PHIS3*(.5*X3*X3-X3*sin(X3)+.5*sin(X3)*...
			   sin(X3))/(W3^4);
Q3(2,1)=Q3(1,2);
Q3(1,3)=PHIS3*(sin(X3)-X3*cos(X3)-.5*X3+.25*...
			   sin(2*X3))/(W3^3);
Q3(3,1)=Q3(1,3);
Q3(1,4)=PHIS3*(cos(X3)+X3*sin(X3)-.5*sin(X3)*...
			   sin(X3)-1)/(W3*W3);
Q3(4,1)=Q3(1,4);
Q3(2,2)=PHIS3*(1.5*X3-2*sin(X3)+.25*sin(2*X3))/(W3^3);
Q3(2,3)=PHIS3*(1-cos(X3)-.5*sin(X3)*sin(X3))/(W3*W3);
Q3(3,2)=Q3(2,3);
Q3(2,4)=PHIS3*(sin(X3)-.5*X3-.25*sin(2*X3))/W3;
Q3(4,2)=Q3(2,4);
Q3(3,3)=PHIS3*(.5*X3-.25*sin(2*X3))/W3;
Q3(3,4)=.5*PHIS3*sin(X3)*sin(X3);
Q3(4,3)=Q3(3,4);
Q3(4,4)=W3*PHIS3*(.5*X3+.25*sin(2*X3));

P1(1,1)=SIGN2;
P1(2,2)=(VM*HEDEGFIL/57.3)^2;
P1(3,3)=XNT*XNT;
P1(4,4)=W1*W1*XNT*XNT;
	
P2(1,1)=SIGN2;
P2(2,2)=(VM*HEDEGFIL/57.3)^2;
P2(3,3)=XNT*XNT;
P2(4,4)=W2*W2*XNT*XNT;
	
P3(1,1)=SIGN2;
P3(2,2)=(VM*HEDEGFIL/57.3)^2;
P3(3,3)=XNT*XNT;
P3(4,4)=W3*W3*XNT*XNT;
	
HMAT=[1 0 0 0];
HT=HMAT';
PHIT1=PHI1';
PHIT2=PHI2';
PHIT3=PHI3';
T=0.;
H=.001;
S=0.;
XNC=0.;
XNL=0.;
XLAM=Y/RTM;
YTDD=XNTREAL*sin(WREAL*T);
YTDDD=XNTREAL*WREAL*cos(WREAL*T);
	
YH1=0.;
YDH1=0.;
YTDDH1=0.;
YTDDDH1=0.;
	
YH2=0.;
YDH2=0.;
YTDDH2=0.;
YTDDDH2=0.;
	
YH3=0.;
YDH3=0.;
YTDDH3=0.;
YTDDDH3=0.;
	
PROB1=.333;
PROB2=.333;
PROB3=.333;
while T<=(TF-.0001)
	S=S+H;
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
		TGO=TF-T+.000001;
		RTM=VC*TGO;
		XLAM=Y/(VC*TGO);
		YTDD=XNTREAL*sin(WREAL*T);
		XNLD=(XNC-XNL)/TAU;
		YDD=YTDD-XNL;
		FLAG=1;
	end
	FLAG=0;
 	Y=.5*(YOLD+Y+H*YD);
 	YD=.5*(YDOLD+YD+H*YDD);
	XNL=.5*(XNLOLD+XNL+H*XNLD);
	if S>=(TS-.00001)
 		S=0.;
		TGO=TF-T+.000001;
		RTM=VC*TGO;
		SIGPOS=RTM*SIGNOISE;
		SIGN2=SIGPOS^2;
		PHIP1=PHI1*P1;
		PHIPPHIT1=PHIP1*PHIT1;
		M1=PHIPPHIT1+Q1;
 		HM1=HMAT*M1;
 		HMHT1=HM1*HT;
		HMHTR1=HMHT1(1,1)+SIGN2;
		HMHTRINV1(1,1)=1./HMHTR1;
		MHT1=M1*HT;
		GAIN1=MHT1*HMHTRINV1;
		KH1=GAIN1*HMAT;
		IKH1=IDNP-KH1;
		P1=IKH1*M1;
		
		PHIP2=PHI2*P2;
		PHIPPHIT2=PHIP2*PHIT2;
		M2=PHIPPHIT2+Q2;
 		HM2=HMAT*M2;
 		HMHT2=HM2*HT;
		HMHTR2=HMHT2(1,1)+SIGN2;
		HMHTRINV2(1,1)=1./HMHTR2;
		MHT2=M2*HT;
		GAIN2=MHT2*HMHTRINV2;
		KH2=GAIN2*HMAT;
		IKH2=IDNP-KH2;
		P2=IKH2*M2;
		
		PHIP3=PHI3*P3;
		PHIPPHIT3=PHIP3*PHIT3;
		M3=PHIPPHIT3+Q3;
 		HM3=HMAT*M3;
 		HMHT3=HM3*HT;
		HMHTR3=HMHT3(1,1)+SIGN2;
		HMHTRINV3(1,1)=1./HMHTR3;
		MHT3=M3*HT;
		GAIN3=MHT3*HMHTRINV3;
		KH3=GAIN3*HMAT;
		IKH3=IDNP-KH3;
		P3=IKH3*M3;
 	
		CPZ1=HMHTR1;
		CPZ2=HMHTR2;
		CPZ3=HMHTR3;
	
		YTDD=XNTREAL*sin(WREAL*T);
		YTDDD=XNTREAL*WREAL*cos(WREAL*T);
		XLAMNOISE=SIGNOISE*randn;
		YSTAR=RTM*(XLAM+XLAMNOISE);
	
		RES1=YSTAR-YH1-TS*YDH1-(1-cos(X1))*YTDDH1/...
			(W1*W1)-(X1-sin(X1))*YTDDDH1/(W1*W1*W1)...
			+.5*TS*TS*XNL;
		YH1=YH1+TS*YDH1+(1-cos(X1))*YTDDH1/(W1*W1)+...
			(X1-sin(X1))*YTDDDH1/(W1*W1*W1)+GAIN1(1,1)*RES1-...
			.5*TS*TS*XNL;
		YDH1=YDH1+sin(X1)*YTDDH1/W1+(1-cos(X1))*YTDDDH1/(W1*W1)...
			+GAIN1(2,1)*RES1-TS*XNL;
     		YTDDHNEW1=cos(X1)*YTDDH1+sin(X1)*YTDDDH1/W1+...
			GAIN1(3,1)*RES1;
		YTDDDH1=-W1*sin(X1)*YTDDH1+cos(X1)*YTDDDH1+GAIN1(4,1)*RES1;
		YTDDH1=YTDDHNEW1;
	
		RES2=YSTAR-YH2-TS*YDH2-(1-cos(X2))*YTDDH2/(W2*W2)-...
			(X2-sin(X2))*YTDDDH2/(W2*W2*W2)+.5*TS*TS*XNL;
		YH2=YH2+TS*YDH2+(1-cos(X2))*YTDDH2/(W2*W2)+...
			(X2-sin(X2))*YTDDDH2/(W2*W2*W2)+...
			GAIN2(1,1)*RES2-.5*TS*TS*XNL;
		YDH2=YDH2+sin(X2)*YTDDH2/W2+(1-cos(X2))*...
			YTDDDH2/(W2*W2)+GAIN2(2,1)*RES2-TS*XNL;
		YTDDHNEW2=cos(X2)*YTDDH2+sin(X2)*YTDDDH2/W2+...
			GAIN2(3,1)*RES2;
		YTDDDH2=-W2*sin(X2)*YTDDH2+cos(X2)*YTDDDH2...
			+GAIN2(4,1)*RES2;
		YTDDH2=YTDDHNEW2;
	
		RES3=YSTAR-YH3-TS*YDH3-(1-cos(X3))*YTDDH3/...
			(W3*W3)-(X3-sin(X3))*YTDDDH3/(W3*W3*W3)...
			+.5*TS*TS*XNL;
		YH3=YH3+TS*YDH3+(1-cos(X3))*YTDDH3/(W3*W3)+...
			(X3-sin(X3))*YTDDDH3/(W3*W3*W3)+...
			GAIN3(1,1)*RES3-.5*TS*TS*XNL;
		YDH3=YDH3+sin(X3)*YTDDH3/W3+(1-cos(X3))*...
			YTDDDH3/(W3*W3)+GAIN3(2,1)*RES3-TS*XNL;
		YTDDHNEW3=cos(X3)*YTDDH3+sin(X3)*YTDDDH3/W3+...
			GAIN3(3,1)*RES3;
		YTDDDH3=-W3*sin(X3)*YTDDH3+cos(X3)*YTDDDH3...
			+GAIN3(4,1)*RES3;
		YTDDH3=YTDDHNEW3;
	
		F1=exp(-.5*RES1*RES1/CPZ1)/sqrt(6.28*CPZ1);
		F2=exp(-.5*RES2*RES2/CPZ2)/sqrt(6.28*CPZ2);
		F3=exp(-.5*RES3*RES3/CPZ3)/sqrt(6.28*CPZ3);
	
		PROB1=PROB1*F1/(PROB1*F1+PROB2*F2+PROB3*F3);
		PROB2=PROB2*F2/(PROB1*F1+PROB2*F2+PROB3*F3);
		PROB3=PROB3*F3/(PROB1*F1+PROB2*F2+PROB3*F3);
	
		WHPZ=W1*PROB1+W2*PROB2+W3*PROB3;
		YHPZ=YH1*PROB1+YH2*PROB2+YH3*PROB3;
		YDHPZ=YDH1*PROB1+YDH2*PROB2+YDH3*PROB3;
		YTDDHPZ=YTDDH1*PROB1+YTDDH2*PROB2+YTDDH3*PROB3;
		YTDDDHPZ=YTDDDH1*PROB1+YTDDDH2*PROB2+...
			YTDDDH3*PROB3;
	
		XS=TGO/TAU;
		TOP=6.*XS*XS*(exp(-XS)-1.+XS);
		BOT1=2*XS*XS*XS+3.+6.*XS-6.*XS*XS;
		BOT2=-12.*XS*exp(-XS)-3.*exp(-2.*XS);
		XNPP=TOP/(.0001+BOT1+BOT2);
		C1=XNPP/(TGO*TGO);
		C2=XNPP/TGO;
		C3=XNPP*(1.-cos(WHPZ*TGO))/(WHPZ*WHPZ*TGO*TGO);
		C4=-XNPP*(exp(-XS)+XS-1.)/(XS*XS);
		C5=XNPP*(WHPZ*TGO-sin(WHPZ*TGO))/(WHPZ*...
				WHPZ*WHPZ*TGO*TGO);
		XNC=C1*YHPZ+C2*YDHPZ+C3*YTDDHPZ+C4*XNL+...
				C5*YTDDDHPZ;
		if XNC>XLIM
			XNC=XLIM;
		end
		if XNC<-XLIM
			XNC=-XLIM;
		end
		YTDDG=YTDD/32.2;
		YTDDHPZG=YTDDHPZ/32.2;
		count=count+1;
		ArrayT(count)=T;
		ArrayPROB1(count)=PROB1;
		ArrayPROB2(count)=PROB2;
		ArrayPROB3(count)=PROB3;
		ArrayWREAL(count)=WREAL;
		ArrayWHPZ(count)=WHPZ;
		ArrayYTDDG(count)=YTDDG;
		ArrayYTDDHPZG(count)=YTDDHPZG;
	end
end
output=[ArrayT',ArrayPROB1',ArrayPROB2',ArrayPROB3',...
		ArrayWREAL',ArrayWHPZ',ArrayYTDDG',...
		ArrayYTDDHPZG'];
save datfil.txt output /ascii
disp 'simulation finished'
clc
figure
plot(ArrayT,ArrayPROB1,ArrayT,ArrayPROB2,ArrayT,...
		ArrayPROB3),grid
xlabel('Time (s) ')
ylabel('Probability')	
axis([0 10 0 1.2])
figure
plot(ArrayT,ArrayWREAL,ArrayT,ArrayWHPZ),grid
xlabel('Time (s) ')
ylabel('Frequency (r/s)')	
axis([0 10 0 4])
figure
plot(ArrayT,ArrayYTDDG,ArrayT,ArrayYTDDHPZG),grid
xlabel('Time (s) ')
ylabel('Acceleration (g)')	
axis([0 10 -6 6])
