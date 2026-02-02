clear
XNT=96.6;
XNP=4.;
TAU=1.;
TF=10.;
VM=3000.;
HEDEG=-20.;
APN=0;
T=0.;
S=0.;
TP=T+.00001;
X1=0.;
X2=0.;
X3=1.;
X4=0.;
XNPP=0.;
H=.01;
HE=HEDEG/57.3;
n=0.;
while TP<=(TF-1e-5)
	X1OLD=X1;
	X2OLD=X2;
	X3OLD=X3;
	X4OLD=X4;
   	STEP=1;
   	FLAG=0;
	while STEP<=1
      		if FLAG==1
         		STEP=2;
         		X1=X1+H*X1D;
         		X2=X2+H*X2D;
         		X3=X3+H*X3D;
         		X4=X4+H*X4D;
         		TP=TP+H;
      		end
      		TGO=TP+.00001;
      		if APN==0
         		C1=XNP/(TGO*TGO);
         		C2=XNP/TGO;
         		C3=0.;
         		C4=0.;
      		elseif APN==1
         		C1=XNP/(TGO*TGO);
         		C2=XNP/TGO;
         		C3=.5*XNP;
         		C4=0.;
      		else
         		X=TGO/TAU;
         		TOP=6.*X*X*(exp(-X)-1.+X);
         		BOT1=2*X*X*X+3.+6.*X-6.*X*X;
         		BOT2=-12.*X*exp(-X)-3.*exp(-2.*X);
         		XNPP=TOP/(.0001+BOT1+BOT2);
         		C1=XNPP/(TGO*TGO);
         		C2=XNPP/TGO;
         		C3=.5*XNPP;
         		C4=-XNPP*(exp(-X)+X-1.)/(X*X);
      		end
      		X1D=X2+C3*X4/TAU;
      		X2D=X3+C2*X4/TAU;
      		X3D=C1*X4/TAU;
      		X4D=-X4/TAU-X2+C4*X4/TAU;
      		FLAG=1;
   	end
   	FLAG=0;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	X2=(X2OLD+X2)/2+.5*H*X2D;
	X3=(X3OLD+X3)/2+.5*H*X3D;
   	X4=(X4OLD+X4)/2+.5*H*X4D;
   	S=S+H;
	if S>=.0999
      		S=0.;
      		n=n+1;
      		ArrayTP(n)=TP;
      		ArrayXMNT(n)=XNT*X1;
      		ArrayXMHE(n)=-VM*HE*X2;
   	end
end
figure
plot(ArrayTP,ArrayXMNT),grid
xlabel('Flight Time (Sec)')
ylabel('Target Maneuver Miss (Ft)')
clc
output=[ArrayTP',ArrayXMNT',ArrayXMHE'];
save datfil.txt output  -ascii
disp 'simulation finished'


