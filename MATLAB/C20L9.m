clear
n=0;
VC=4000.;
XNT=0.;
DISPLACE=200.;
VM=3000.;
TAU=1.;
XNP=3.;
XNCLIM=99999999.;
TF=10.;
for THOM = 0.1:0.1:10
	QSWITCH=0;
	Y=0.;
	YD=0.;
	XNL=0.;
	D=0.;
	ELAMDH=0.;
	X4=0.;
	X5=0.;
	T=0.;
	H=.01;
	S=0.;
	while ~(T>(TF-.0001))
		TGO=TF-T+.00001;
		if TGO<=THOM & QSWITCH==0
			QSWITCH=1;
			Y=Y+DISPLACE;
			XLAM=Y/(VC*TGO);
			D=XLAM;
		end
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
			TGO=TF-T+.00001;
			XLAM=Y/(VC*TGO);
			DD=5.*(XLAM-D)/TAU;
			ELAMDHD=5.*(DD-ELAMDH)/TAU;
			XNC=XNP*VC*ELAMDH;
			if XNC>XNCLIM
				XNC=XNCLIM;
			end
			if XNC<-XNCLIM
				XNC=-XNCLIM;
			end
			X4D=5.*(XNC-X4)/TAU;
			X5D=5.*(X4-X5)/TAU;
			XNLD=5.*(X5-XNL)/TAU;
			YDD=XNT-XNL;
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
	ArrayTHOM(n)=THOM;
	ArrayY(n)=Y;
end
figure
plot(ArrayTHOM,ArrayY),grid
xlabel('Homing Time (s)')
ylabel('Miss (Ft)')
clc
output=[ArrayTHOM',ArrayY'];
save datfil.txt output  -ascii
disp 'simulation finished'