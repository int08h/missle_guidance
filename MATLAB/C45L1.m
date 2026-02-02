clear
count=0;
IOPTION=0;
ITGT=1;
XNTAV=117.6;
if ITGT==1
	TF=180.;
else
	TF=240.;
end
PRED=0.*3280.;
VM=9000.;
VC=18000.;
XNTMAX=9.*32.2;
XNCMAX=966.;
APN=0.;
HEDEG=-57.3*PRED/(VM*TF);
YD=-VM*HEDEG/57.3;
Y=0.;
XNP=3.;
T=0.;
H=.001;
S=0.;
DELV=0.;
SUM=0.;
XN=0.;
while T<(TF-.0001)
	YOLD=Y;
	YDOLD=YD;
	DELVOLD=DELV;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
			STEP=2;
			Y=Y+H*YD;
			YD=YD+H*YDD;
			DELV=DELV+H*DELVD;
			T=T+H;
		end;
		TGO=TF-T+.00001;
		if ITGT==1
			if IOPTION==0
				XNT=XNTMAX*(T/TF)^2;
			else
				if T<180.
					WGT=-212.*T+44000.;
					TRST=54100.;
				else
					WGT=3300.;
					TRST=0.;
				end
					XNT=32.2*TRST/WGT;
			end
		else
			if IOPTION==0
				XNT=XNTAV;
			else
				if T<120.
					WGT=-2622*T+440660.;
					TRST=725850.;
				elseif T<240.
					WGT=-642.*T+168120.;
					TRST=182250.;
				else
					WGT=5500.;
					TRST=0.;
				end
				XNT=32.2*TRST/WGT;
			end
		end
		XLAMD=(Y+YD*TGO)/(VC*TGO*TGO);
		XNC=XNP*VC*XLAMD+.5*APN*XNP*XNT;
		if XNC>XNCMAX
			XNC=XNCMAX;
		end
		if XNC<-XNCMAX
			XNC=-XNCMAX;
		end
		DELVD=abs(XNC);
		YDD=XNT-XNC;
		FLAG=1;
	end
	FLAG=0;
	Y=.5*(YOLD+Y+H*YD);
	YD=.5*(YDOLD+YD+H*YDD);
	DELV=.5*(DELVOLD+DELV+H*DELVD);
	S=S+H;
	if S>=.09999
		S=0.;
		SUM=SUM+XNT;
		XN=XN+1.;
		count=count+1;
		ArrayT(count)=T;
		ArrayXNT(count)=XNT/32.2;
		ArrayXNC(count)=XNC/32.2;
		ArrayDELV(count)=DELV/3.28;
	end
end
figure
plot(ArrayT,ArrayXNC),grid 
xlabel('Missile Flight Time (s)')
ylabel('KKV Acceleration (g)  ')
clc
output=[ArrayT',ArrayXNT',ArrayDELV'];
save datfil.txt output -ascii
disp 'simulation finished'
Y
DELV/3.28
