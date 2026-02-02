clear
count=0;
VC=4000.;
TAP=.2;
XNT=96.6;
BETA=.8;
XNP=3.;
TS=.1;
TS2=.02;
for TF=.1:.1:10.0,
	Y=0.;
	YD=0.;
	T=0.;
	H=.001;
	S=0.;
	S2=0.;
	GFILTER=1.-BETA^2;
	HFILTER=(1.-BETA)^2;
	XLAMHOLD=0.;
	XLAMDHOLD=0.;
	Y1OLD=0.;
	Y2OLD=0.;
	Y3OLD=0.;
	Y4OLD=0.;
	Y5OLD=0.;
	XNC=0.;
	XNL=0.;
	Y1NEW=0.;
	PZ=0.;
	XLAM=0.;
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
			XNLD=(XNC-XNL)/TAP;
			YDD=XNT-XNL;
			FLAG=1;
		end;
		FLAG=0;
 		Y=.5*(YOLD+Y+H*YD);
 		YD=.5*(YDOLD+YD+H*YDD);
		XNL=.5*(XNLOLD+XNL+H*XNLD);
		S=S+H;
		S2=S2+H;
		if S2>=(TS2 - 1e-5)
			S2=0.;
			Y1NEW=XLAM;
			Y2NEW=Y1OLD;
			Y3NEW=Y2OLD;
			Y4NEW=Y3OLD;
			Y5NEW=Y4OLD;
			PZ=.2*(Y5OLD+Y5NEW+Y4NEW+Y3NEW+Y2NEW+XLAM);
			Y1OLD=Y1NEW;
			Y2OLD=Y2NEW;
			Y3OLD=Y3NEW;
			Y4OLD=Y4NEW;
			Y5OLD=Y5NEW;
		end;
		if S>=(TS - 1e-5)
 			S=0.;
			RES=PZ-(XLAMHOLD+TS*XLAMDHOLD);
			XLAMHNEW=GFILTER*RES+XLAMHOLD+TS*XLAMDHOLD;
			XLAMDHNEW=HFILTER*RES/TS+XLAMDHOLD;
			XNC=XNP*VC*XLAMDHNEW;
			XLAMHOLD=XLAMHNEW;
			XLAMDHOLD=XLAMDHNEW;
		end;
	end;
	count=count+1;
	ArrayTF(count)=TF;
	ArrayY(count)=Y;
end;
figure
plot(ArrayTF',ArrayY'),grid
title('Standard miss for various flight times')
xlabel('Flight Time (S)')
ylabel('Miss (Ft) ')
axis([00,10,-40,100])
clc
output=[ArrayTF',ArrayY'];
save datfil.txt output -ascii
disp('Simulation Complete')
