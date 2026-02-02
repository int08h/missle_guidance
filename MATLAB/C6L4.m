clear
XNP=3.;
TAU=1.;
TF=10.;
VC=1.;
PHIFN=1;
PHIRN=1;
PHIRNA=1;
PHIGL=1;
RA=1;
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
X9=0.;
X10=0.;
X11=0.;
X12=0.;
H=.01;
n=0.;
while TP<=(TF-1e-5)
      X2OLD=X2;
      X3OLD=X3;
      X4OLD=X4;
      X5OLD=X5;
      X6OLD=X6;
      X7OLD=X7;
      X8OLD=X8;
      X9OLD=X9;
      X10OLD=X10;
      X11OLD=X11;
      X12OLD=X12;
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
        X9=X9+H*X9D;
        X10=X10+H*X10D;
                        X11=X11+H*X11D;
        X12=X12+H*X12D;
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
        X9D=Y1^2;
        X10D=(Y1*VC*TGO/RA)^2;
                X11D=(Y1*((VC*TGO/RA)^2))^2;
                X12D=(Y1/(VC*TGO))^2;
        FLAG=1;
      end
      FLAG=0;
      X2=.5*(X2OLD+X2+H*X2D);
      X3=.5*(X3OLD+X3+H*X3D);
      X4=.5*(X4OLD+X4+H*X4D);
      X5=.5*(X5OLD+X5+H*X5D);
      X6=.5*(X6OLD+X6+H*X6D);
      X7=.5*(X7OLD+X7+H*X7D);
      X8=.5*(X8OLD+X8+H*X8D);
      X9=.5*(X9OLD+X9+H*X9D);
      X10=.5*(X10OLD+X10+H*X10D);
      X11=.5*(X11OLD+X11+H*X11D);
      X12=.5*(X12OLD+X12+H*X12D);
      S=S+H;
      if S>=.0999
        S=0.;
        n=n+1;
        XMFN=sqrt(X9*PHIFN);
        XMRN=sqrt(X10*PHIRN);
        XMRNA=sqrt(X11*PHIRNA);
        XMGL=sqrt(X12*PHIGL);
        ArrayTP(n)=TP;
        ArrayXMFN(n)=XMFN;
        ArrayXMRN(n)=XMRN;
        ArrayXMRNA(n)=XMRNA;
                ArrayXMGL(n)=XMGL;
      end
end
XMFN
XMRN
XMRNA
XMGL
figure
plot(ArrayTP,ArrayXMFN),grid
xlabel('Normalized Flight Time (Sec)')
ylabel('Normalized Fading Noise Miss')
figure
plot(ArrayTP,ArrayXMRN),grid
xlabel('Normalized Flight Time (Sec)')
ylabel('Missile Semiactive Noise Miss')
figure
plot(ArrayTP,ArrayXMRNA),grid
xlabel('Normalized Flight Time (Sec)')
ylabel('Missile Active Noise Miss')
figure
plot(ArrayTP,ArrayXMGL),grid
xlabel('Normalized Flight Time (Sec)')
ylabel('Missile Glint Noise Miss')
clc
output=[ArrayTP',ArrayXMFN',ArrayXMRN',ArrayXMRNA',ArrayXMGL'];
save datfil.txt output -ascii
disp 'simulation finished'
