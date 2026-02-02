clear all
n=0;
XNP=3.;
TAU=1.;
TF=10.;
YTIC=1.;
VC=4000.;
T=0.;
S=0.;
TP=T+.00001;
X2=0;
X3=1;
X4=0;
X5=0.;
X6=0.;
X7=0.;
X8=0.;
H=.01;
while ~(TP>(TF-.00001))
	S=S+H;
	X2OLD=X2;
	X3OLD=X3;
	X4OLD=X4;
	X5OLD=X5;
	X6OLD=X6;
	X7OLD=X7;
	X8OLD=X8;
	STEP=1;
	FLAG=0;
	while STEP<=1
		if FLAG==1
			STEP=2;
			X2=X2+H*X2D;
			X3=X3+H*X3D;
			X4=X4+H*X4D;
			X5=X5+H*X5D;
			X6=X6+H*X6D;
			X7=X7+H*X7D;
			X8=X8+H*X8D;
			TP=TP+H;
		end
		X2D=X3;
		Y1=5.*(5.*X5/TAU+X4)/TAU;
		TGO=TP+.00001;
		X3D=Y1/(VC*TGO);
		X4D=-Y1;
		X5D=-5.*X5/TAU+5.*X6*XNP*VC/TAU;
		X6D=-5.*X6/TAU+5.*X7/TAU;
		X7D=-5.*X7/TAU+5.*X8/TAU;
		X8D=-5.*X8/TAU-X2;
		FLAG=1;
	end
	FLAG=0;
	X2=(X2OLD+X2)/2+.5*H*X2D;
	X3=(X3OLD+X3)/2+.5*H*X3D;
	X4=(X4OLD+X4)/2+.5*H*X4D;
	X5=(X5OLD+X5)/2+.5*H*X5D;
	X6=(X6OLD+X6)/2+.5*H*X6D;
	X7=(X7OLD+X7)/2+.5*H*X7D;
	X8=(X8OLD+X8)/2+.5*H*X8D;
	if S<.09999
		S=0.;
		XMYT=YTIC*X3;
		n=n+1;
		ArrayTP(n)=TP;
		ArrayXMYT(n)=XMYT;
	end
end
plot(ArrayTP,ArrayXMYT),grid
xlabel('Homing Time (s)')
ylabel('Miss (Ft)')
clc
output=[ArrayTP',ArrayXMYT'];
save datfil.txt output  -ascii
disp 'simulation finished'
