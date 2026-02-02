clear
count=0;
XNT = 96.6;
XNP = 3;
TAU = 1;
TF = 10;
T=0;
S=0;
TP=T+.00001;
X1=0;
X2=0;
X3=1;
X4=0;
X5=0.;
X6=0;
X7=0;
H=.01;
XNU=.5;
BETA=96.6
while TP <= (TF - 1e-5)
	STEP=1;
	FLAG=0;
	S=S+H;
	X1OLD=X1;
	X2OLD=X2;
	X3OLD=X3;
	X4OLD=X4;
	X5OLD=X5;
	X6OLD=X6;
	X7OLD=X7;
	while STEP <=1
		if FLAG==1
			STEP=2;
			X1=X1+H*X1D;
			X2=X2+H*X2D;
			X3=X3+H*X3D;
			X4=X4+H*X4D;
			X5=X5+H*X5D;
			X6=X6+H*X6D;
			X7=X7+H*X7D;
			TP=TP+H;
		end;
		X1D=X2;
		X2D=X3;
		Y1=(X4-X2)/TAU;
		TGO=TP+.00001;
		X3D=XNP*Y1/TGO;
		X4D=-Y1;
		X5D=X1*X1;
		X6D=X2-2.*XNU*X6; 
		X7D=(2.*XNU*X6)^2;
		FLAG=1;
	end;
	FLAG=0;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	X2=(X2OLD+X2)/2+.5*H*X2D;
	X3=(X3OLD+X3)/2+.5*H*X3D;
	X4=(X4OLD+X4)/2+.5*H*X4D;
   	X5=(X5OLD+X5)/2+.5*H*X5D;
	X6=(X6OLD+X6)/2+.5*H*X6D;
	X7=(X7OLD+X7)/2+.5*H*X7D;
   	S=S+H;
   	if S>=.000999
      		S=0.;
			XMBT=BETA*sqrt(X7/XNU);
			XIC=BETA*X6;
			RMS=sqrt(XMBT^2+XIC^2);
      		count=count+1;
      		ArrayTP(count)=TP;
      		ArrayRMS(count)=RMS;
   	end;
end
figure
plot(ArrayTP, ArrayRMS),grid
title('Adjoint model using shaping filter approach')
xlabel('Flight Time (S)')
ylabel('RMS Miss (Ft)')
%axis([00,10,00,30])
clc
output=[ArrayTP',ArrayRMS'];
save datfil.txt output -ascii
disp('Simulation Complete')

