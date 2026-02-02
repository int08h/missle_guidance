clear
VC=4000.;
XNT=0.;
Y=0.;
VM=3000.;
HEDEG=-20.;
TF=10.;
XNP=4.;
YD=-VM*HEDEG/57.3;
T=0.;
H=.01;
S=0.;
n=0.;
while T<=(TF-1e-5)
   YOLD=Y;
   YDOLD=YD;   
   STEP=1;
   FLAG=0;
   while STEP<=1
      if FLAG==1
         STEP=2;
         Y=Y+H*YD;
         YD=YD+H*YDD;
         T=T+H;
      end
      TGO=TF-T+.00001;
      XLAMD=(Y+YD*TGO)/(VC*TGO*TGO);
      XNC=XNP*VC*XLAMD;
      YDD=XNT-XNC;
      FLAG=1;
   end
   FLAG=0;
   Y=.5*(YOLD+Y+H*YD);
   YD=.5*(YDOLD+YD+H*YDD);
   S=S+H;
   if S>=.0999
      S=0.;
      n=n+1;
      ArrayT(n)=T;
      ArrayY(n)=Y;
      ArrayYD(n)=YD;
      ArrayXNCG(n)=XNC/32.2;
   end
end
figure
plot(ArrayT,ArrayXNCG),grid
xlabel('Time (Sec)')
ylabel('Missile Acceleration (G)')
clc
output=[ArrayT',ArrayY',ArrayYD',ArrayXNCG'];
save datfil.txt output  -ascii
disp 'simulation finished'
