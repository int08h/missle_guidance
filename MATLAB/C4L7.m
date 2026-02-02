clear
n=0;
XNT=193.2;
XNP=3.;
TAU=.5;
TF=10.;
VM=3000.;
HEDEG=-20.;
VC=4000.;
TS=.1;
RA=30000.;
SIGGL=10.;
SIGFN=.002;
SIGRN=.02;
SIGRNA=.02;
PHIGL=SIGGL*SIGGL*TS;
PHIFN=SIGFN*SIGFN*TS;
PHIRN=SIGRN*SIGRN*TS;
PHIRNA=SIGRNA*SIGRNA*TS;
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
	S=S+H;
   	if S>=.0999
      	S=0.;
      	n=n+1;
		XMGL=sqrt(PHIGL*X5);
		XMFN=sqrt(PHIFN*X6);
		XMRN=sqrt(PHIRN*X7);
		XMRNA=sqrt(PHIRNA*X8);
		ArrayTP(n)=TP;
		ArrayXMGL(n)=XMGL;
        ArrayXMGLTH(n)=1.44*sqrt(PHIGL/TAU);
		ArrayXMFN(n)=XMFN;
        ArrayXMFNTH(n)=.532*VC*sqrt(TAU*PHIFN);
		ArrayXMRN(n)=XMRN;
        ArrayXMRNTH(n)=1.06*VC*VC*(TAU^1.5)*sqrt(PHIRN)/RA;
		ArrayXMRNA(n)=XMRNA;
        ArrayXMRNATH(n)=4.66*VC*VC*VC*(TAU^2.5)*sqrt(PHIRNA)/(RA*RA);
	end
end
figure
plot(ArrayTP,ArrayXMGL,ArrayTP,ArrayXMGLTH),grid
xlabel('Flight Time (Sec)')
ylabel('Glint Miss Standard Deviation (Ft)')
figure
plot(ArrayTP,ArrayXMFN,ArrayTP,ArrayXMFNTH),grid
xlabel('Flight Time (Sec)')
ylabel('Range Independent Noise Miss Stndard Deviation(Ft)')
figure
plot(ArrayTP,ArrayXMRN,ArrayTP,ArrayXMRNTH),grid
xlabel('Flight Time (Sec)')
ylabel('Semiactive Range Dependent Noise Miss Stndard Deviation(Ft)')
figure
plot(ArrayTP,ArrayXMRNA,ArrayTP,ArrayXMRNATH),grid
xlabel('Flight Time (Sec)')
ylabel('Active Range Dependent Noise Miss Stndard Deviation(Ft)')
clc
output=[ArrayTP',ArrayXMGL',ArrayXMGLTH',ArrayXMFN',ArrayXMFNTH',...
    ArrayXMRN',ArrayXMRNTH',ArrayXMRNA',ArrayXMRNATH'];
save datfil output  -ascii
disp 'simulation finished'

