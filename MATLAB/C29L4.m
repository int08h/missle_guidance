clear
n=0;
VC=4000.;
XNT=193.2;
XNP=3.;
XNCLIM=99999999.;
TAU=.25;
W=2.;
WH=2.;
APN=1;
for TF=.1:.1:10
	Y=0.;
	YD=0.;
	XNL=0.;
	D=0.;
	ELAMDH=0.;
	X4=0.;
	X5=0.;
	T=0.;
	H=.01;
	while T<=(TF-.0001)
 		YOLD=Y;
		YDOLD=YD;
		XNLOLD=XNL;
		DOLD=D;
		ELAMDHOLD=ELAMDH;
		X4OLD=X4;
		X5OLD=X5;
		STEP=1;
		FLAG=0;
      		while STEP<=1
         		if FLAG==1
            			STEP=2;
 				Y=Y+H*YD;
 				YD=YD+H*YDD;
				XNL=XNL+H*XNLD;
				ELAMDH=ELAMDH+H*ELAMDHD;
				D=D+H*DD;
				X4=X4+H*X4D;
				X5=X5+H*X5D;
				T=T+H;
			end
			YTDD=XNT*sin(W*T);
			YTDDD=W*XNT*cos(W*T);
 			TGO=TF-T+.00001;
			XLAM=Y/(VC*TGO);
			DD=5.*(XLAM-D)/TAU;
			ELAMDHD=5.*(DD-ELAMDH)/TAU;
			if APN==1
				XNC=XNP*VC*ELAMDH;
			elseif APN==2
				XP=WH*TGO;
            			XNC=XNP*VC*ELAMDH+XNP*YTDD*(1.-cos(XP))/XP^2+...
               			XNP*YTDDD*(XP-sin(XP))/(XP*XP*WH);
			else
				X=TGO/TAU;
				XP=WH*TGO;
				TOP=6.*X*X*(exp(-X)-1.+X);
				BOT1=2*X*X*X+3.+6.*X-6.*X*X;
				BOT2=-12.*X*exp(-X)-3.*exp(-2.*X);
				XNPP=TOP/(.0001+BOT1+BOT2);
            			XNC=XNPP*VC*ELAMDH+XNPP*YTDD*(1.-cos(XP))/XP^2....
               				+XNPP*YTDDD*(XP-sin(XP))/(XP*XP*WH)-...
               			XNPP*XNL*TAU*TAU*(exp(-X)+X-1.)/TGO^2;
			end
			if XNC>XNCLIM
				XNC=XNCLIM;
			end
			if XNC<-XNCLIM
				XNC=-XNCLIM;
			end
			X4D=5.*(XNC-X4)/TAU;
			X5D=5.*(X4-X5)/TAU;
			XNLD=5.*(X5-XNL)/TAU;
			YDD=YTDD-XNL;
			FLAG=1;
      		end
      		FLAG=0;
 		Y=.5*(YOLD+Y+H*YD);
 		YD=.5*(YDOLD+YD+H*YDD);
		XNL=.5*(XNLOLD+XNL+H*XNLD);
		D=.5*(DOLD+D+H*DD);
		ELAMDH=.5*(ELAMDHOLD+ELAMDH+H*ELAMDHD);
		X4=.5*(X4OLD+X4+H*X4D);
		X5=.5*(X5OLD+X5+H*X5D);
	end
	n=n+1;
	ArrayTF(n)=TF;
	ArrayY(n)=Y;
end
figure
plot(ArrayTF,ArrayY),grid
xlabel('Flight Time (Sec)')
ylabel('Miss (Ft)')
clc
output=[ArrayTF',ArrayY'];
save datfil.txt output  -ascii
disp 'simulation finished'
