clear
n=0;	
TS=.1;
APN=0;
XLIMG=5.;
VM=3000.;
BETA=1000.;
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
VTM1=0.-VM1;
VTM2=0.-VM2;
VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
XLIM=XLIMG*32.2;
XNC=0.;
while ~(T>0.&RM2<=0.)
	if RTM<1000.
		H=.0001;
	else
		H=.01;
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
		VTM1=0.-VM1;
		VTM2=0.-VM2;
		VC=-(RTM1*VTM1+RTM2*VTM2)/RTM;
		XLAM=atan2(RTM2,RTM1);
		XLAMD=(RTM1*VTM2-RTM2*VTM1)/(RTM*RTM);
 		if RM2<30000.
			RHO=.002378*exp(-RM2/30000);
		else
			RHO=.0034*exp(-RM2/22000);
		end
		VM=sqrt(VM1^2+VM2^2);
		Q=.5*RHO*VM*VM;
		GAM=atan2(VM2,VM1);
		DRAG=Q*32.2/BETA;
		XNE1=-XNC*sin(XLAM);
		XNE2=XNC*cos(XLAM);
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
	if S>=(TS-.00001)
		S=0.;
		RM1K=RM1/1000.;
		RM2K=RM2/1000.;
		RT1K=RT1/1000.;
		RT2K=RT2/1000.;
		DRAG1=-DRAG*cos(GAM);
		DRAG2=-DRAG*sin(GAM)-32.2;
		DRAGPLOS=-DRAG1*sin(XLAM)+DRAG2*cos(XLAM);
		ATPLOS=0.;
		if T>30.
			if APN==0
				XNC=3.*VC*XLAMD;
			else
				XNC=3.*VC*XLAMD+1.5*(ATPLOS-DRAGPLOS);
			end
		else
			XNC=0.;
		end
		if XNC>XLIM
			XNC=XLIM;
		end
		if XNC<-XLIM
			XNC=-XLIM;
		end
		XNCG=XNC/32.2;
		DRAGPLOSG=DRAGPLOS/32.2;
		n=n+1;
		ArrayT(n)=T;
		ArrayRT1K(n)=RT1K;
		ArrayRT2K(n)=RT2K;
		ArrayRM1K(n)=RM1K;
		ArrayRM2K(n)=RM2K;
		ArrayXNCG(n)=XNC/32.2;
	end
end
RTM
figure
plot(ArrayRT1K,ArrayRT2K,ArrayRM1K,ArrayRM2K),grid
xlabel('Downrange (Ft) ')
ylabel('Altitude (Ft)')
figure
plot(ArrayT,ArrayXNCG),grid
xlabel('Time (sec)')
ylabel('Acceleration of missle (G)')
clc
output=[ArrayT',ArrayRT1K',ArrayRT2K',ArrayRM1K',ArrayRM2K',ArrayXNCG'];
save datfil output -ascii
disp '*** Simulation Complete'
	
	
