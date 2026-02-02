clear
n=0;
TAU=.5;
APN=0;
ORDER=4;
MVR=1;
VC=9000.;
W=2.;
WREAL=2.;
WH=W;
XNT=96.6;
XNTREAL=96.6;
TS=.01;
YIC=0.;
VM=3000.;
HEDEG=0.;
HEDEGFIL=20.;
XNP=3.;
SIGRIN=.001;
SIGGL=0.;
RA=21000.;
SRN=0.;
TF=10.;
QPERFECT=0;
PHASE=0./57.3;
X=WH*TS;
Y=YIC;
YD=-VM*HEDEG/57.3;
PHIS=WH*WH*XNT*XNT/TF;
RTM=VC*TF;
SIGNOISE=sqrt(SIGRIN^2+(SIGGL/RTM)^2+(SRN*RTM*RTM/(RA*RA))^2);
SIGPOS=RTM*SIGNOISE;
SIGN2=SIGPOS^2;
PHI=zeros(ORDER);
P=zeros(ORDER);
Q=zeros(ORDER);
IDNP=eye(ORDER);
PHI(1,1)=1;
PHI(1,2)=TS;
PHI(1,3)=(1-cos(X))/(WH*WH);
PHI(1,4)=(X-sin(X))/(WH*WH*WH);
PHI(2,2)=1;
PHI(2,3)=sin(X)/WH;
PHI(2,4)=(1-cos(X))/(WH*WH);
PHI(3,3)=cos(X);
PHI(3,4)=sin(X)/WH;
PHI(4,3)=-WH*sin(X);
PHI(4,4)=cos(X);
Q(1,1)=PHIS*(.333*X^3-2*sin(X)+2*X*cos(X)+.5*X-.25*sin(2*X))/(WH^5);
Q(1,2)=PHIS*(.5*X*X-X*sin(X)+.5*sin(X)*sin(X))/(WH^4);
Q(2,1)=Q(1,2);
Q(1,3)=PHIS*(sin(X)-X*cos(X)-.5*X+.25*sin(2*X))/(WH^3);
Q(3,1)=Q(1,3);
Q(1,4)=PHIS*(cos(X)+X*sin(X)-.5*sin(X)*sin(X)-1)/(WH*WH);
Q(4,1)=Q(1,4);
Q(2,2)=PHIS*(1.5*X-2*sin(X)+.25*sin(2*X))/(WH^3);
Q(2,3)=PHIS*(1-cos(X)-.5*sin(X)*sin(X))/(WH*WH);
Q(3,2)=Q(2,3);
Q(2,4)=PHIS*(sin(X)-.5*X-.25*sin(2*X))/WH;
Q(4,2)=Q(2,4);
Q(3,3)=PHIS*(.5*X-.25*sin(2*X))/WH;
Q(3,4)=.5*PHIS*sin(X)*sin(X);
Q(4,3)=Q(3,4);
Q(4,4)=WH*PHIS*(.5*X+.25*sin(2*X));
P(1,1)=SIGN2;
P(2,2)=(VM*HEDEGFIL/57.3)^2;
P(3,3)=XNT*XNT;
P(4,4)=WH*WH*XNT*XNT;
HMAT=[1 0 0 0];
HT=HMAT';
PHIT=PHI';	
T=0.;
H=.001;
S=0.;
XNC=0.;
XNL=0.;
XLAM=Y/RTM;
if MVR==0
	YTDD=XNTREAL;
	YTDDD=0.;
else
	YTDD=XNTREAL*sin(WREAL*T);
	YTDDD=XNTREAL*WREAL*cos(WREAL*T);
end
if QPERFECT==1
	YH=Y;
	YDH=YD;
	YTDDH=YTDD;
	YTDDDH=YTDDD;
else
	YH=0.;
	YDH=0.;
	YTDDH=0.;
	YTDDDH=0.;
end
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
		TGO=TF-T+.000001;
		RTM=VC*TGO;
		XLAM=Y/(VC*TGO);
		if MVR==0
			YTDD=XNTREAL;
		else
			YTDD=XNTREAL*sin(WREAL*T);
		end
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
		TGO=TF-T+.000001;
		RTM=VC*TGO;
		SIGNOISE=sqrt(SIGRIN^2+(SIGGL/RTM)^2+(SRN*RTM*RTM/(RA*RA))^2);
		SIGPOS=RTM*SIGNOISE;
		SIGN2=SIGPOS^2;
		RMAT=[SIGN2];
		PHIP=PHI*P;
		PHIPPHIT=PHIP*PHIT;
		M=PHIPPHIT+Q;
		HM=HMAT*M;
		HMHT=HM*HT;
 		HMHTR=HMHT+RMAT;
		HMHTRINV=inv(HMHTR);
		MHT=M*HT;
		GAIN=MHT*HMHTRINV;
		KH=GAIN*HMAT;
		IKH=IDNP-KH;
		P=IKH*M;
 		if MVR==0
			YTDD=XNTREAL;
			YTDDD=0.;
		else
			YTDD=XNTREAL*sin(WREAL*T);
			YTDDD=XNTREAL*WREAL*cos(WREAL*T);
		end
		XLAMNOISE=SIGNOISE*randn;
		YSTAR=RTM*(XLAM+XLAMNOISE);
		RES=YSTAR-YH-TS*YDH-(1-cos(X))*YTDDH/(WH*WH)-(X-sin(X))...
				*YTDDDH/(WH*WH*WH)+.5*TS*TS*XNL;
		YH=YH+TS*YDH+(1-cos(X))*YTDDH/(WH*WH)+(X-sin(X))...
				*YTDDDH/(WH*WH*WH)+GAIN(1,1)*RES-.5*TS*TS*XNL;
		YDH=YDH+sin(X)*YTDDH/WH+(1-cos(X))*YTDDDH/(WH*WH)...
				+GAIN(2,1)*RES-TS*XNL;
		YTDDHNEW=cos(X)*YTDDH+sin(X)*YTDDDH/WH+GAIN(3,1)*RES;
		YTDDDH=-WH*sin(X)*YTDDH+cos(X)*YTDDDH+GAIN(4,1)*RES;
		YTDDH=YTDDHNEW;
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
		YTDDG=YTDD/32.2	;
		YTDDHG=YTDDH/32.2;	
		ERRY=Y-YH;
		SP11=sqrt(P(1,1));
		SP11P=-SP11;
		ERRYD=YD-YDH;
		SP22=sqrt(P(2,2));
		SP22P=-SP22;
		ERRYTDDG=(YTDD-YTDDH)/32.2;
		SP33G=sqrt(P(3,3))/32.2;
		SP33GN=-SP33G;
		ERRYTDDDG=(YTDDD-YTDDDH)/32.2;
		SP44G=sqrt(P(4,4))/32.2;
		SP44GN=-SP44G;
		YTDDG=YTDD/32.2;
		YTDDHG=YTDDH/32.2;
		YTDDDG=YTDDD/32.2;
		YTDDDHG=YTDDDH/32.2;
		n=n+1;
		ArrayT(n)=T;
		ArrayYTDDG(n)=YTDDG;
		ArrayYTDDHG(n)=YTDDHG;
		ArrayYTDDDG(n)=YTDDDG;
		ArrayYTDDDHG(n)=YTDDDHG;
		ArrayERRYTDDG(n)=ERRYTDDG;
		ArraySP33G(n)=SP33G;
		ArraySP33GN(n)=SP33GN;
		ArrayERRYTDDDG(n)=ERRYTDDDG;
		ArraySP44G(n)=SP44G;
		ArraySP44GN(n)=SP44GN;
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
plot(ArrayT,ArrayERRYTDDG,ArrayT,ArraySP33G,ArrayT,ArraySP33GN)...
	,grid
xlabel('Time (Sec)')
ylabel('Error in Estimate of Acceleration (G)')
figure
plot(ArrayT,ArrayERRYTDDDG,ArrayT,ArraySP44G,ArrayT,ArraySP44GN)...
	,grid
xlabel('Time (Sec)')
ylabel('Error in Estimate of Jerk (G/S)')
clc
output=[ArrayT',ArrayYTDDG',ArrayYTDDHG',ArrayYTDDDG',...
		ArrayYTDDDHG'];
save datfil.txt output  -ascii
output=[ArrayT',ArrayERRYTDDG',ArraySP33G',ArraySP33GN'...
		,ArrayERRYTDDDG',ArraySP44G',ArraySP44GN'];
save covfil.txt output  -ascii
disp 'simulation finished'	
