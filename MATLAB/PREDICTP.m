function[RM1F]= PREDICTP(TP,RM1P,RM2P,VM1P,VM2P,XNC1F,RT1P,RT2P,BETAH)
H=.01;
RM1=RM1P;
RM2=RM2P;
VM1=VM1P;
VM2=VM2P;
XNC=XNC1F;
RT1=RT1P;
RT2=RT2P;
BETA=BETAH;
T=TP;
RTM1=RT1-RM1;
RTM2=RT2-RM2;
RTM=sqrt(RTM1^2+RTM2^2);
while RM2>=0.
	if RTM<1000.
%Integration interval increased by factor of 10 to increase speed
		%H=.0001
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
end
RM1F=RM1;
