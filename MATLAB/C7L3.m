clear
count=1;
XNT=96.6;
XNP=3.;
TF=10.;
TS=.1;
BETA=.8;
SIGNOISE=.001;
VC=4000.;
T=0.;
S=0.;
TP=T+.00001;
X1=0;
X2=0;
X3=1;
X5=0.;
Y1OLD=0.;
Y2OLD=0.;
Y3OLD=0.;
Y4OLD=0.;
Y5OLD=0.;
H=.01;
GFILTER=1.-BETA^2;
HFILTER=(1.-BETA)^2;
while TP<=(TF-1e-5)	
   	X1OLD=X1;
	X2OLD=X2;
	X3OLD=X3;
	X5OLD=X5;
	STEP=1;
	FLAG=0;
   	while STEP<=1
      		if FLAG==1
         		STEP=2;
         		X1=X1+H*X1D;
         		X2=X2+H*X2D;
         		X3=X3+H*X3D;
         		X5=X5+H*X5D;
         		TP=TP+H;
      		end
      		TGO=TP+.00001;
      		X1D=X2;
      		X2D=X3+Y4OLD/(VC*TGO);
      		X3D=(Y4OLD)/(VC*TGO*TGO);
      		X5D=-X2;
      		FLAG=1;
   	end
   	FLAG=0;
	X1=(X1OLD+X1)/2+.5*H*X1D;
	X2=(X2OLD+X2)/2+.5*H*X2D;
	X3=(X3OLD+X3)/2+.5*H*X3D;
   	X5=(X5OLD+X5)/2+.5*H*X5D;
   	S=S+H;
   	if S>=(TS-.0001)
      		S=0.;
      		TEMP1=(X5-Y1OLD)*XNP*VC;
      		TEMP2=HFILTER*(Y2OLD+TEMP1)/TS+GFILTER*Y3OLD;
      		Y1NEW=X5;
      		Y2NEW=TEMP1+Y2OLD+TS*(Y3OLD-TEMP2);
      		Y3NEW=Y3OLD-TEMP2;
      		Y4NEW=Y4OLD+TEMP2;
      		Y5NEW=Y5OLD+TEMP2*TEMP2;
      		Y1OLD=Y1NEW;
      		Y2OLD=Y2NEW;
      		Y3OLD=Y3NEW;
      		Y4OLD=Y4NEW;
      		Y5OLD=Y5NEW;
      		XMNOISE=SIGNOISE*sqrt(Y5NEW);
      		XMNT=XNT*X1;
      		count=count+1;
      		ArrayTP(count)=TP;
      		ArrayXMNT(count)=XMNT;
      		ArrayXMNOISE(count)=XMNOISE;
   	end
end
figure
plot(ArrayTP,ArrayXMNT),grid
xlabel('Flight Time (Sec)')
ylabel('Target Maneuver Miss (Ft)')
figure
plot(ArrayTP,ArrayXMNOISE),grid
xlabel('Flight Time (Sec)')
ylabel('Noise Miss (Ft)')
clc
output=[ArrayTP',ArrayXMNT',ArrayXMNOISE'];
save datfil.txt output  -ascii
disp 'simulation finished'



