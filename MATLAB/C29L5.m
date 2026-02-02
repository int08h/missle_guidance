clear all
close all

XNT=96.6;
XNP=4;
TAU=1.;
VM=3000;
HEDEG=-20;

PI=3.1416;
TF=25.6;
TS=.1;
XNPOINT=TF/TS;
NPOINT=round(XNPOINT);
FS=1./TS;
T=0.;
S=0.;
TP=T+.00001;
X1=0;
X2=0;
X3=1;
X4=0;
H=.01;
HE=HEDEG/57.3;
J=0;
while TP<=(TF-.0001)
	X1OLD=X1;
	X2OLD=X2;
	X3OLD=X3;
	X4OLD=X4;
	STEP=1;
	FLAG=0;
	while STEP <=1
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
	if S>=(TS-.00001)
		S=0;
        J=J+1;
		XMNT=XNT*X1;
		XMHE=-VM*HE*X2;
		X(J)=X2; 
    end
end



for K=1:NPOINT
    XR(K)=0;
    XI(K)=0;
    for N=1:NPOINT
        AG    = 2*PI*K*N/NPOINT;
        XR(K) = XR(K)+X(N)*cos(AG)/NPOINT;
        XI(K) = XI(K)-X(N)*sin(AG)/NPOINT;
    end
end

IMAX  = NPOINT/2;
for I=1:IMAX+1
    F=FS*I/NPOINT;
    XMAG(I)=sqrt(XR(I)^2+XI(I)^2);
    W=2.*PI*F;
    PZ=XMAG(I)*XNT*XNPOINT/FS;
    ArrayW(I)=W;
    ArrayPZ(I)=PZ;
end

figure
plot(ArrayW,ArrayPZ),grid
xlabel('W (r/s)')
ylabel('PZ ')
clc
output=[ArrayW',ArrayPZ'];
save datfil.txt output  -ascii
disp 'simulation finished'