clear
n=0;
XNT=193.2;
XNP=3.;
TAU=1.;
TF=10.;
VC=4000.;
W=3.;
T=0.;
S=0.;
TP=T+.00001;
X2=0;
X3=1;
X4=0.;
X5=0.;
X6=0.;
H=.01;
while TP<=(TF-1e-5)
	X2OLD=X2;
	X3OLD=X3;
	X4OLD=X4;
	X5OLD=X5;
	X6OLD=X6;
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
			TP=TP+H;
		end
		X2D=X3;
		Y1=(-X2+X4)/TAU;
		TGO=TP+.00001;
		X3D=Y1*XNP/TGO;
		X4D=-Y1;
		X5D=X2-W*W*X6;
		X6D=X5;
		FLAG=1;
   	end
   	FLAG=0;
   	X2=(X2OLD+X2)/2+.5*H*X2D;
   	X3=(X3OLD+X3)/2+.5*H*X3D;
	X4=(X4OLD+X4)/2+.5*H*X4D;
	X5=(X5OLD+X5)/2+.5*H*X5D;
	X6=(X6OLD+X6)/2+.5*H*X6D;
	S=S+H;
	if S>=.09999
      		S=0.;
      		n=n+1;
		XMWEAVE=XNT*W*X6;
  	   	ArrayTP(n)=TP;
		ArrayXMWEAVE(n)=XMWEAVE;
     	end
end
figure
plot(ArrayTP,ArrayXMWEAVE),grid
xlabel('Flight Time (S)')
ylabel('Miss (Ft) ')
clc
output=[ArrayTP',ArrayXMWEAVE'];
save datfil.txt output /ascii
disp 'simulation finished'
