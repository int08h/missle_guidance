n=0;
XNT=193.2;
XNP=3.;
TAU=.5;
TF=10.;
VM=3000.;
HEDEG=-20.;
VC=4000.;
TS=.1;
SIGRN=.01;
RA=30000.;
SIGGL=2.;
SIGFN=.001;
SIGRN=.01;
SIGRNA=.01;
PHIGL=SIGGL*SIGGL/TS;
PHIFN=SIGFN*SIGFN/TS;
PHIRN=SIGRN*SIGRN/TS;
PHIRNA=SIGRNA*SIGRNA/TS;
RA=30000.;
T=0.;
S=0.;
TP=T+.00001;
X1=0;
X2=0;
X3=1;
X4=0;
X5=0.;
X6=0.;
X7=0.;
X8=0.;
X9=0.;
H=.01;
HE=HEDEG/57.3;
while TP<=(TF-1e-5)
	S=S+H;
	X1OLD=X1;
	X2OLD=X2;
	X3OLD=X3;
	X4OLD=X4;
	X5OLD=X5;
	X6OLD=X6;
	X7OLD=X7;
	X8OLD=X8;
	X9OLD=X9;
    STEP=1;
	FLAG=0;
	while STEP<=1
      		if FLAG==1
         	STEP=2;
			X1=X1+H*X1D;
			X2=X2+H*X2D;
			X3=X3+H*X3D;
			X4=X4+H*X4D;
			X5=X5+H*X5D;
			X6=X6+H*X6D;
			X7=X7+H*X7D;
			X8=X8+H*X8D;
			X9=X9+H*X9D;
			TP=TP+H;
            end
        X1D=X2;
		X2D=X3;
		Y1=(X4-X2)/TAU;
		TGO=TP+.00001;
		X3D=XNP*Y1/TGO;
		X4D=-Y1;
		X5D=X3D^2;
		X6D=(XNP*VC*Y1)^2;
		X7D=(XNP*VC*VC*Y1*TGO/RA)^2;
		X8D=(XNP*VC*VC*VC*Y1*TGO*TGO/RA^2)^2;
		X9D=X1^2;
		FLAG=1;
   	end
   	FLAG=0;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	X2=(X2OLD+X2)/2+.5*H*X2D;
	X3=(X3OLD+X3)/2+.5*H*X3D;
	X4=(X4OLD+X4)/2+.5*H*X4D;
	X5=(X5OLD+X5)/2+.5*H*X5D;
	X6=(X6OLD+X6)/2+.5*H*X6D;
	X7=(X7OLD+X7)/2+.5*H*X7D;
	X8=(X8OLD+X8)/2+.5*H*X8D;
	X9=(X9OLD+X9)/2+.5*H*X9D;
	S=S+H;
   	if S>=.0999
      	S=0.;
      	n=n+1;
		XMGL=sqrt(PHIGL*X5);
		XMFN=sqrt(PHIFN*X6);
		XMRN=sqrt(PHIRN*X7);
		XMRNA=sqrt(PHIRNA*X8);
		XMUDNT=XNT*sqrt(X9/TGO);
		RMSSA=sqrt(XMGL^2+XMFN^2+XMRN^2+XMUDNT^2);
		RMSA=sqrt(XMGL^2+XMFN^2+XMRNA^2+XMUDNT^2);
		ArrayTP(n)=TP;
		ArrayXMGL(n)=XMGL;
		ArrayXMFN(n)=XMFN;
		ArrayXMRN(n)=XMRN;
		ArrayXMRNA(n)=XMRNA;
		ArrayXMUDNT(n)= XMUDNT;
		ArrayRMSSA(n)= RMSSA;
		ArrayRMSA(n)=RMSA;
	end
end
figure
plot(ArrayTP,ArrayXMGL,ArrayTP,ArrayXMFN,ArrayTP,ArrayXMRN,ArrayTP,...
    ArrayXMUDNT,ArrayTP,ArrayRMSSA),grid
xlabel('Flight Time (Sec)')
ylabel('RMS Miss For Semiactive Error Budget (Ft)')
figure
plot(ArrayTP,ArrayXMGL,ArrayTP,ArrayXMFN,ArrayTP,ArrayXMRNA,ArrayTP,...
    ArrayXMUDNT,ArrayTP,ArrayRMSA),grid
xlabel('Flight Time (Sec)')
ylabel('RMS Miss For Active Error Budget (Ft)')
clc
output=[ArrayTP',ArrayXMGL',ArrayXMFN',ArrayXMRN',ArrayXMRNA',...
    ArrayRMSSA',ArrayRMSA'];
save datfil output  -ascii
disp 'simulation finished'

