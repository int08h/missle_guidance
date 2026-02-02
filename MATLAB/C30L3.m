clear
n=0;
PHIS2=0.;
XNT=96.6;
W=2.;
PHASEDEG=0.;
SIGRIN=.0001;
SIGGL=0.;
SRN=0.;
RA=21000.;
WHIC=-1.;
TS=.01;
TF=10.;
PHIS1=W*W*XNT*XNT/TF;
QPERFECT=0;
VC=9000.;
XNP=3.;
XNCLIM=9999999.;
APN=4;
TAU=.5;
HEDEG=0.;
VM=3000.;
QEKF=0;
PHASE=PHASEDEG/57.3;
ORDER=5;
TGO=TF;
T=0.;
X=W*T;
S=0.;
Y=0.;
YD=-XNT/W-VM*HEDEG/57.3;
YTDD=XNT*sin(W*T);
YTDDD=XNT*W*cos(W*T);
XNC=0.;
XNL=0.;
H=.001;
HP=.001;
TS2=TS*TS;
TS3=TS2*TS;
TS4=TS3*TS;
TS5=TS4*TS;
TS6=TS5*TS;
TS7=TS6*TS;
WH=WHIC;
if QPERFECT==1
	YH=Y;
	YDH=YD;
	YTDDH=YTDD;
	YTDDDH=YTDDD;
	WH=W;
else
	YH=0.;
	YDH=0.;
	YTDDH=0.;
	YTDDDH=0.;
end
PHI=zeros(ORDER);
P=zeros(ORDER);
Q=zeros(ORDER);
IDNP=eye(ORDER);
RTM=VC*TF;
SIGNOISE=sqrt(SIGRIN^2+(SIGGL/RTM)^2+(SRN*RTM*RTM/(RA*RA))^2);
YNOISE=SIGNOISE*RTM;
P(1,1)=YNOISE*YNOISE;
P(2,2)=(VM*20./57.3)^2;
P(3,3)=XNT*XNT;
P(4,4)=(W*XNT)^2;
P(5,5)=W^2;
HMAT=[1 0 0 0 0];
HT=HMAT';
while T<=(TF-.0001)
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
		YTDD=XNT*sin(W*T);
 		TGO=TF-T+.00001;
		XNLD=(XNC-XNL)/TAU;
		YDD=YTDD-XNL;
		FLAG=1;
	end
	FLAG=0;
 	Y=.5*(YOLD+Y+H*YD);
 	YD=.5*(YDOLD+YD+H*YDD);
	XNL=.5*(XNLOLD+XNL+H*XNLD);
	S=S+H;
	if S>=(TS-.00001)
		S=0.;
		YTDD=XNT*sin(W*T);
		YTDDD=XNT*W*cos(W*T);
		PHI(1,1)=1.;
		PHI(1,2)=TS;
		PHI(2,2)=1.;
		PHI(2,3)=TS;
		PHI(3,3)=1.;
		PHI(3,4)=TS;
		PHI(4,3)=-WH*WH*TS;
		PHI(4,4)=1.;
		PHI(4,5)=-2.*WH*YTDDH*TS;
		PHI(5,5)=1.;
		Q(3,3)=PHIS1*TS*TS*TS/3.;
		Q(3,4)=PHIS1*TS*TS/2.;
		Q(4,3)=Q(3,4);
		Q(4,4)=4.*WH*WH*YTDDH*YTDDH*PHIS2*TS*TS*TS/3.+PHIS1*TS;
		Q(4,5)=-WH*YTDDH*TS*TS*PHIS2;
		Q(5,4)=Q(4,5);
		Q(5,5)=PHIS2*TS;
		PHIT=PHI';
		PHIP=PHI*P;
		PHIPPHIT=PHIP*PHIT;
		M=PHIPPHIT+Q;
		HM=HMAT*M;
 		HMHT=HM*HT;
 		RTM=VC*TGO;
		SIGNOISE=sqrt(SIGRIN^2+(SIGGL/RTM)^2+(SRN*RTM*RTM/...
			(RA*RA))^2);
		YNOISE=SIGNOISE*RTM;
		RMAT=[YNOISE^2];
		HMHTR=HMHT+RMAT;
		HMHTRINV=inv(HMHTR);
		MHT=M*HT;
		GAIN=MHT*HMHTRINV;
		KH=GAIN*HMAT;
		IKH=IDNP-KH;
 		P=IKH*M;
 		RTM=VC*TGO;
		XLAM=Y/RTM;
		XNOISE=SIGNOISE*randn;
		XLAMS=XLAM+XNOISE;
		[YB,YDB,YTDDB,YTDDDB]=PROJECT(T,TS,YH,YDH,YTDDH,YTDDDH,...
			HP,XNL,WH);
		RES=RTM*XLAMS-YB;
		YH=YB+GAIN(1,1)*RES;
		YDH=YDB+GAIN(2,1)*RES;
		YTDDH=YTDDB+GAIN(3,1)*RES;
		YTDDDH=YTDDDB+GAIN(4,1)*RES;
		WH=WH+GAIN(5,1)*RES;
		if APN==0
			XNC=XNP*(YH+YDH*TGO)/(TGO*TGO);
		elseif APN==1
			XNC=XNP*(YH+YDH*TGO+.5*YTDDH*TGO*TGO)/(TGO*TGO);
		elseif APN==2
			XS=TGO/TAU;
			TOP=6.*XS*XS*(exp(-XS)-1.+XS);
			BOT1=2*XS*XS*XS+3.+6.*XS-6.*XS*XS;
			BOT2=-12.*XS*exp(-XS)-3.*exp(-2.*XS);
			XNPP=TOP/(.0001+BOT1+BOT2);
			C1=XNPP/(TGO*TGO);
			C2=XNPP/TGO;
			C3=.5*XNPP;
			C4=-XNPP*(exp(-XS)+XS-1.)/(XS*XS);
			XNC=C1*YH+C2*YDH+C3*YTDDH+C4*XNL;
		elseif APN==3
			XP=WH*TGO;
			XNC=XNP*(YH+YDH*TGO)/(TGO*TGO)+XNP*YTDDH*...
			(1.-cos(XP))/XP^2+XNP*YTDDDH*(XP-sin(XP))/(XP*XP*WH);
		else
			XS=TGO/TAU;
			TOP=6.*XS*XS*(exp(-XS)-1.+XS);
			BOT1=2*XS*XS*XS+3.+6.*XS-6.*XS*XS;
			BOT2=-12.*XS*exp(-XS)-3.*exp(-2.*XS);
			XNPP=TOP/(.0001+BOT1+BOT2);
			C1=XNPP/(TGO*TGO);
			C2=XNPP/TGO;
			C3=XNPP*(1.-cos(WH*TGO))/(WH*WH*TGO*TGO);
			C4=-XNPP*(exp(-XS)+XS-1.)/(XS*XS);
			C5=XNPP*(WH*TGO-sin(WH*TGO))/(WH*WH*WH*TGO*TGO);
			XNC=C1*YH+C2*YDH+C3*YTDDH+C4*XNL+C5*YTDDDH;
		end
		if XNC>XNCLIM
			XNC=XNCLIM;
		end
		if XNC<-XNCLIM
			XNC=-XNCLIM;
		end
		ERRYTDD=YTDD-YTDDH;
		ERRYTDDG=ERRYTDD/32.2;
		ERRYTDDD=YTDDD-YTDDDH;
		ERRYTDDDG=ERRYTDDD/32.2;
		ERRW=W-WH;
		SP44=sqrt(P(4,4));
		SP44P=-SP44;
		SP33=sqrt(P(3,3));
		SP33P=-SP33;
		SP33G=SP33/32.2;
		SP33PG=SP33P/32.2;
		SP44G=SP44/32.2;
		SP44PG=SP44P/32.2;
		SP55=sqrt(P(5,5));
		SP55P=-SP55;
		YTDDG=YTDD/32.2;
		YTDDHG=YTDDH/32.2;
		YTDDDG=YTDDD/32.2;
		YTDDDHG=YTDDDH/32.2;
		XNCG=XNC/32.2;
		n=n+1;
		ArrayT(n)=T;
		ArrayYTDDG(n)=YTDDG;
		ArrayYTDDHG(n)=YTDDHG;
		ArrayYTDDDG(n)=YTDDDG;
		ArrayYTDDDHG(n)=YTDDDHG;
		ArrayW(n)=W;
		ArrayWH(n)=WH;
		ArrayERRYTDDG(n)=ERRYTDDG;
		ArraySP33G(n)=SP33G;
		ArraySP33PG(n)=SP33PG;
		ArrayERRYTDDDG(n)=ERRYTDDDG;
		ArraySP44G(n)=SP44G;
		ArraySP44PG(n)=SP44PG;
		ArrayERRW(n)=ERRW;
		ArraySP55(n)=SP55;
		ArraySP55P(n)=SP55P;
	end
end
figure
plot(ArrayT,ArrayYTDDG,ArrayT,ArrayYTDDHG),grid
xlabel('Time (Sec)')
ylabel('Acceleration and Estimate (G)')
figure
plot(ArrayT,ArrayYTDDDG,ArrayT,ArrayYTDDDHG),grid
xlabel('Time (Sec)')
ylabel('Jerk and Estimate (G/S)')
figure
plot(ArrayT,ArrayW,ArrayT,ArrayWH),grid
xlabel('Time (Sec)')
ylabel('Frequency and Estimate (G/S)')
figure
plot(ArrayT,ArrayERRYTDDG,ArrayT,ArraySP33G,ArrayT,...
	 ArraySP33PG),grid
xlabel('Time (Sec)')
ylabel('Error in Estimate of Acceleration (G)')
figure
plot(ArrayT,ArrayERRYTDDDG,ArrayT,ArraySP44G,ArrayT,...
	 ArraySP44PG),grid
xlabel('Time (Sec)')
ylabel('Error in Estimate of Jerk (G/S)')
clc
output=[ArrayT',ArrayYTDDG',ArrayYTDDHG',ArrayYTDDDG',...
		ArrayYTDDDHG',ArrayW',ArrayWH'];
save datfil.txt output  -ascii
output=[ArrayT',ArrayERRYTDDG',ArraySP33G',ArraySP33PG',...
		ArrayERRYTDDDG',ArraySP44G',ArraySP44PG'];
save covfil.txt output  -ascii
disp 'simulation finished'
