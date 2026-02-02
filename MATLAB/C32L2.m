clear
n=0;
H=.01;
VM=3000.;
BETA=1000.;
BETAH=1000.;
XLIMG=5.;
T=0.;
S=0.;
GAMDEG=45.;
VM1=VM*cos(GAMDEG/57.3);
VM2=VM*sin(GAMDEG/57.3);
RM1=0.;
RM2=0.;
XNC=0.;
RT1=60000.;
RT2=0.;
RTM1=RT1-RM1;
RTM2=RT2-RM2;
RTM=sqrt(RTM1^2+RTM2^2);
XLIM=XLIMG*32.2;
while ~(T>0.&RM2<=0.)
	if RTM<1000.
%Integration interval increased by factor of 10 to increase speed
		%H=.0001;
		H=.001;
	else
		%H=.01
		H=.1;
	end
 	RM1OLD=RM1;
	RM2OLD=RM2;
	VM1OLD=VM1;
	VM2OLD=VM2;
	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
         		STEP=2;
 			RM1=RM1+H*VM1;
			RM2=RM2+H*VM2;
			VM1=VM1+H*AM1;
			VM2=VM2+H*AM2;
			T=T+H;
		end
		RTM1=RT1-RM1;
		RTM2=RT2-RM2;
		RTM=sqrt(RTM1^2+RTM2^2);
		XLAM=atan2(RTM2,RTM1);
		XNE1=-XNC*sin(XLAM);
		XNE2=XNC*cos(XLAM);
 		if RM2<30000.
			RHO=.002378*exp(-RM2/30000);
		else
			RHO=.0034*exp(-RM2/22000);
		end
		VM=sqrt(VM1^2+VM2^2);
		Q=.5*RHO*VM*VM;
		GAM=atan2(VM2,VM1);
		DRAG=Q*32.2/BETA;
 		AM1=-DRAG*cos(GAM)+XNE1;
		AM2=-32.2-DRAG*sin(GAM)+XNE2;
		FLAG=1;
	end
	FLAG=0;
 	RM1=.5*(RM1OLD+RM1+H*VM1);
	RM2=.5*(RM2OLD+RM2+H*VM2);
 	VM1=.5*(VM1OLD+VM1+H*AM1);
	VM2=.5*(VM2OLD+VM2+H*AM2);
	S=S+H;
	if S>=.99999
		S=0.;
		if T>30.
			[X1]=PREDICTP(T,RM1,RM2,VM1,VM2,XNC,RT1,RT2,BETAH);
			DELX=RT1-X1;
			XNCP=XNC+1.;
			[X2]=PREDICTP(T,RM1,RM2,VM1,VM2,XNCP,RT1,RT2,BETAH);
			DXDNC=(X2-X1)/(XNCP-XNC);
			DELXNC=DELX/DXDNC;
			XNC=XNC+DELXNC;
			if XNC>XLIM
				XNC=XLIM;
			end
			if XNC<-XLIM
				XNC=-XLIM;
			end
		end
		RM1K=RM1/1000.;
		RM2K=RM2/1000.;
		RT1K=RT1/1000.;
		RT2K=RT2/1000.;
		XNCG=XNC/32.2;
		n=n+1;
		ArrayT(n)=T;
		ArrayRM1K(n)=RM1K;
		ArrayRM2K(n)=RM2K;
		ArrayRT1K(n)=RT1K;
		ArrayRT2K(n)=RT2K;
		ArrayXNCG(n)=XNCG;
	end
end
RTM
figure
plot(ArrayRM1K,ArrayRM2K,ArrayRT1K,ArrayRT2K),grid
xlabel('Downrange (Kft)')
ylabel('Altitude (KFt)')
figure
plot(ArrayT,ArrayXNCG),grid
xlabel('Time (Sec)')
ylabel('Acceleration (G)')
clc
output=[ArrayT',ArrayRM1K',ArrayRM2K',ArrayRT1K',ArrayRT2K',ArrayXNCG'];
save datfil output  -ascii
disp 'simulation finished'
	
	
