clear
n=0;
XT=-4000.;
YT=5000.;
TF=100.;
H=.01;
TS=.1;
XNC=12.;
XM=0.;
YM=0.;
XMD=0.;
YMD=0.;
S=0.;
T=0.;
PHI=0.;
while T<=TF
	XMOLD=XM;
	YMOLD=YM;
	XMDOLD=XMD;
	YMDOLD=YMD;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
         		STEP=2;
			XM=XM+H*XMD;
			YM=YM+H*YMD;
			XMD=XMD+H*XMDD;
			YMD=YMD+H*YMDD;
			T=T+H;
		end
		XMDD=XNC*cos(PHI);
		YMDD=XNC*sin(PHI);
		FLAG=1;
	end
	FLAG=0;
	XM=(XMOLD+XM)/2+.5*H*XMD;
	YM=(YMOLD+YM)/2+.5*H*YMD;
	XMD=(XMDOLD+XMD)/2+.5*H*XMDD;
	YMD=(YMDOLD+YMD)/2+.5*H*YMDD;
	S=S+H;
	if S>=(TS-.0001)
 		S=0.;
		TGO=TF-T+.0001;
		RTM1=XT-XM;
		RTM2=YT-YM;
		RTM=sqrt(RTM1^2+RTM2^2);
		VTM1=-XMD;
		VTM2=-YMD;
		PHI=atan2(RTM2+VTM2*TGO,RTM1+VTM1*TGO);
		PHIDEG=PHI*57.3;
		n=n+1;
		ArrayT(n)=T;
		ArrayXM(n)=XM;
		ArrayYM(n)=YM;
		ArrayXT(n)=XT;
		ArrayYT(n)=YT;
		ArrayPHIDEG(n)=PHIDEG;
	end
end
RTM
figure
plot(ArrayXM,ArrayYM,ArrayXT,ArrayYT),grid
xlabel('X (Ft)')
ylabel('Y (Ft)')
figure
plot(ArrayT,ArrayPHIDEG),grid
xlabel('Time (Sec)')
ylabel('Roll Rate (Deg/Sec)')
clc
output=[ArrayT',ArrayXM',ArrayYM',ArrayXT',ArrayYT',ArrayPHIDEG'];
save datfil output  -ascii
disp 'simulation finished'
	
