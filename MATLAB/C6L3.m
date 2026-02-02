clear
VC=4000.;
XNT=32.2;
YIC=0.;
VM=3000.;
HEDEG=0.;
TAU=.5;
XNP=3.;
TA=0.;
R=-.01;
n=0.;
for TF=.1:.1:10
   	Y=YIC;
	YD=-VM*HEDEG/57.3;
	YDIC=YD;
	XNL=0.;
   	ELAMDH=0.;	
   	X4=0.;
	X5=0.;
	TH=0.;
   	THH=0.;	
   	T=0.;	
   	H=.01;
      	while T<=(TF-1e-5)
      		YOLD=Y;
      		YDOLD=YD;
      		XNLOLD=XNL;
      		ELAMDHOLD=ELAMDH;
      		X4OLD=X4;
      		X5OLD=X5;
      		THOLD=TH;
      		THHOLD=THH;
      		STEP=1;
      		FLAG=0;
      		while STEP<=1
         		if FLAG==1
            			STEP=2;
            			Y=Y+H*YD;
            			YD=YD+H*YDD;
            			XNL=XNL+H*XNLD;
            			ELAMDH=ELAMDH+H*ELAMDHD;
            			X4=X4+H*X4D;
            			X5=X5+H*X5D;
            			TH=TH+H*THD;
            			THH=THH+H*THHD;
            			T=T+H;
         		end
         		TGO=TF-T+.00001;
         		XLAM=Y/(VC*TGO);
         		EPS=XLAM-TH-THH+R*THH;
         		DD=5.*EPS/TAU;
         		ELAMDHD=5.*(DD-ELAMDH)/TAU;
         		XNC=XNP*VC*ELAMDH;
         		X4D=5.*(XNC-X4)/TAU;
         		X5D=5.*(X4-X5)/TAU;
         		XNLD=5.*(X5-XNL)/TAU;
         		THD=XNL/VM+TA*XNLD/VM;
         		THHD=DD-THD;
         		YDD=XNT-XNL;    
         		FLAG=1;
      		end
      		FLAG=0;
      		Y=.5*(YOLD+Y+H*YD);
      		YD=.5*(YDOLD+YD+H*YDD);
      		XNL=.5*(XNLOLD+XNL+H*XNLD);
      		ELAMDH=.5*(ELAMDHOLD+ELAMDH+H*ELAMDHD);
      		X4=.5*(X4OLD+X4+H*X4D);
      		X5=.5*(X5OLD+X5+H*X5D);
      		TH=.5*(THOLD+TH+H*THD);
      		THH=.5*(THHOLD+THH+H*THHD);
   	end
   	n=n+1;
   	ArrayTF(n)=TF;
   	ArrayY(n)=Y;
end
figure
plot(ArrayTF,ArrayY),grid
xlabel('Flight Time (Sec)')
ylabel('Miss (Ft)')
clc
output=[ArrayTF',ArrayY'];
save datfil.txt output  -ascii
disp 'simulation finished'
