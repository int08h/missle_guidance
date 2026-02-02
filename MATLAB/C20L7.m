clear
n=0;
XNP=3.;
TAU=1.;
TF=5.;
VC=4000.;
DISP=200.;
VM=3000.;
HE=-DISP/VM;
T=0.;
S=0.;
TP=T+.00001;
X1=0;
X2=0;
X3=1;
X4=0;
H=.01;
while ~(TP>(TF-.00001))
	S=S+H;
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
		Y1=(X4-XNP*VC*X2)/TAU;
		TGO=TP+.00001;
		X3D=Y1/(VC*TGO);
		X4D=-Y1;
		FLAG=1;
	end
	FLAG=0;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	X2=(X2OLD+X2)/2+.5*H*X2D;
	X3=(X3OLD+X3)/2+.5*H*X3D;
	X4=(X4OLD+X4)/2+.5*H*X4D;
	if S<.09999
		S=0.;
		XMY=DISP*X3;
		XIC=X4*DISP/(VC*TGO);
		TOT=XMY+XIC;
		XMHE=-VM*HE*X2;
		XMHEPZ=XMHE/TP;
		n=n+1;
		ArrayTP(n)=TP;
		ArrayTOT(n)=TOT;
		ArrayXMHEPZ(n)=XMHEPZ;
	end	
end
figure
plot(ArrayTP,ArrayTOT,ArrayTP,ArrayXMHEPZ),grid
xlabel('Homing Time (s)')
ylabel('Miss (Ft)')
clc
output=[ArrayTP',ArrayTOT',ArrayXMHEPZ'];
save datfil.txt output  -ascii
disp 'simulation finished'
