clear
count=0;
XNT=96.6;
XNP=3.;
TAU=1.;
TF=10.;
T=0.;
S=0.;
TINT=.5;
MISS=0;
TP=T+.00001+TINT;
X1=0;
X2=0;
X5=0.;
if MISS==1
	X3=1.;
	X4=0.;
else
	X3=XNP/(TAU*TINT);
	X4=-1./TAU;
end
H=.01;
while TP<=(TF - 1e-5)
	X1OLD=X1;
	X2OLD=X2;
	X3OLD=X3;
	X4OLD=X4;
	X5OLD=X5;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
			X1=X1+H*X1D;
			X2=X2+H*X2D;
			X3=X3+H*X3D;
			X4=X4+H*X4D;
			X5=X5+H*X5D;
			TP=TP+H;
			STEP=2;
		end
		X1D=X2;
		X2D=X3;
		Y1=(X4-X2)/TAU;
		TGO=TP+.00001;
		X3D=XNP*Y1/TGO;
		X4D=-Y1;
		X5D=X1*X1;
		FLAG=1;
	end;
	FLAG=0;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	X2=(X2OLD+X2)/2+.5*H*X2D;
	X3=(X3OLD+X3)/2+.5*H*X3D;
	X4=(X4OLD+X4)/2+.5*H*X4D;
	X5=(X5OLD+X5)/2+.5*H*X5D;
	S=S+H;
	if S>=.099999
		S=0.;
		XMUDNT=XNT*sqrt(X5/TGO);
		if MISS==0
			XMUDNT=XMUDNT/32.2;
		end
		count=count+1;
		ArrayTP(count)=TP;
		ArrayXMUDNT(count)=XMUDNT;
	end
end
%figure
plot(ArrayTP,ArrayXMUDNT),grid
xlabel('Flight Time (Sec)')
ylabel('Acceleration (G) ')
clc
output=[ArrayTP',ArrayXMUDNT'];
save datfil.txt output -ascii
disp 'simulation finished'
