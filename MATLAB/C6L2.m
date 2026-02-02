clear
XNT=32.2;
XNP=3.;
TAU=1.;
TF=10.;
VC=4000.;
XNTD=32.2;
XNTDD=32.2;
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
X10=0.;
H=.01;
n=0.;
while TP<=(TF-1e-5)
   X1OLD=X1;
   X2OLD=X2;
   X3OLD=X3;
   X4OLD=X4;
   X5OLD=X5;
   X6OLD=X6;
   X7OLD=X7;
   X8OLD=X8;
   X9OLD=X9;
   X10OLD=X10;
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
         X10=X10+H*X10D;
         TP=TP+H;
      end
      X1D=X2;
      X2D=X3;
      Y1=5.*(5.*X5/TAU+X4)/TAU;
      TGO=TP+.00001;
      X3D=Y1/(VC*TGO);
      X4D=-Y1;
      X5D=-5.*X5/TAU+5.*X6*XNP*VC/TAU;
      X6D=-5.*X6/TAU+5.*X7/TAU;
      X7D=-5.*X7/TAU+5.*X8/TAU;
      X8D=-5.*X8/TAU-X2;
      X9D=X1;
      X10D=X9;
      FLAG=1;
   end
   FLAG=0;
   X1=.5*(X1OLD+X1+H*X1D);
   X2=.5*(X2OLD+X2+H*X2D);
   X3=.5*(X3OLD+X3+H*X3D);
   X4=.5*(X4OLD+X4+H*X4D);
   X5=.5*(X5OLD+X5+H*X5D);
   X6=.5*(X6OLD+X6+H*X6D);
   X7=.5*(X7OLD+X7+H*X7D);
   X8=.5*(X8OLD+X8+H*X8D);
   X9=.5*(X9OLD+X9+H*X9D);
   X10=.5*(X10OLD+X10+H*X10D);
   S=S+H;
   if S>=.0999
      S=0.;
      n=n+1;
      ArrayTP(n)=TP;
      ArrayXMNT(n)=XNT*X1;
      ArrayXMNTD(n)=XNTD*X9;
      ArrayXMNTDD(n)=XNTDD*X10;
   end
end
figure
plot(ArrayTP,ArrayXMNT),grid
xlabel('Normalized Flight Time (Sec)')
ylabel('Missile Miss Due To Step Maneuver(Ft/G-Sec^2)')
figure
plot(ArrayTP,ArrayXMNTD),grid
xlabel('Normalized Flight Time (Sec)')
ylabel('Missile Miss Due To Ramp Maneuver(Ft/G-Sec^3/S)')
figure
plot(ArrayTP,ArrayXMNTDD),grid
xlabel('Normalized Flight Time (Sec)')
ylabel('Missile Miss Due To Parabolic Maneuver(Ft/G-Sec^4/S^2)')
clc
output=[ArrayTP',ArrayXMNT',ArrayXMNTD',ArrayXMNTD'];
save datfil.txt output  -ascii
disp 'simulation finished'
