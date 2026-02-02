clear
n=0;
XNT=96.6;
XNP=4.;
TAU=1.;
TF=10;
VM=3000.;
HEDEG=-20.;
T=0.;
S=0.;
TP=T+.00001;
X1=0;
X2=0;
X3=1;
X4=0;
H=.01;
HE=HEDEG/57.3;
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
      		X1D=X2;
      		X2D=X3;
      		Y1=(X4-X2)/TAU;
      		TGO=TP+.00001;
      		X3D=XNP*Y1/TGO;
      		X4D=-Y1;
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
figure
plot(ArrayTP,ArrayXMHE),grid
xlabel('Flight Time (Sec)')
ylabel('Heading Error Miss (Ft)')
clc
output=[ArrayTP',ArrayXMNT',ArrayXMHE'];
save datfil output  -ascii
disp 'simulation finished'


