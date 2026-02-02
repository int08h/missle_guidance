clear
count=0;
H=.01;
VM=3000.;
BETA=1000.;
T=0.;
S=0.;
GAMDEG=45.;
VM1=VM*cos(GAMDEG/57.3);
VM2=VM*sin(GAMDEG/57.3);
RM1=0.;
RM2=0.;
while ~(T>0. & RM2<=0.)
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
      		if RM2<30000.
         		RHO=.002378*exp(-RM2/30000);
      		else
         		RHO=.0034*exp(-RM2/22000);
      		end
      		VM=sqrt(VM1^2+VM2^2);
      		Q=.5*RHO*VM*VM;
      		GAM=atan2(VM2,VM1);
      		DRAG=Q*32.2/BETA;
      		AM1=-DRAG*cos(GAM);
      		AM2=-32.2-DRAG*sin(GAM);
      		FLAG=1;
   	end
   	FLAG=0;
   	RM1=.5*(RM1OLD+RM1+H*VM1);
   	RM2=.5*(RM2OLD+RM2+H*VM2);
   	VM1=.5*(VM1OLD+VM1+H*AM1);
   	VM2=.5*(VM2OLD+VM2+H*AM2);
   	S=S+H;
   	if S>=.099999
      		S=0.;
      		RM1K=RM1/1000.;
      		RM2K=RM2/1000.;
      		count=count+1;
      		ArrayT(count)=T;
      		ArrayRM1K(count)=RM1K;
      		ArrayRM2K(count)=RM2K;
   	end
end
figure
plot(ArrayRM1K,ArrayRM2K),grid
xlabel('Downrange (Kft)')
ylabel('Altitude (Kft) ')
clc
output=[ArrayT',ArrayRM1K',ArrayRM2K'];
save datfil.txt output -ascii
disp 'simulation finished'
