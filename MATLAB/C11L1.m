clear
count=0;
RT1=0.;
RT2=100000.;
VT=6000.;
GAMTDEG=45.;
BETA=500.;
VT1=VT*cos(GAMTDEG/57.3);
VT2=-VT*sin(GAMTDEG/57.3);
T=0.;
H=.01;
S=0.;
while RT2>=0.
	RT1OLD=RT1;
	RT2OLD=RT2;
	VT1OLD=VT1;
	VT2OLD=VT2;
 	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG==1
			RT1=RT1+H*VT1;
			RT2=RT2+H*VT2;
			VT1=VT1+H*AT1;
			VT2=VT2+H*AT2;
 			T=T+H;
 			STEP=2;
 		end
 		if RT2<=30000.
			RHO=.002378*exp(-RT2/30000.);
		else
			RHO=.0034*exp(-RT2/22000.);
		end
		VT=sqrt(VT1^2+VT2^2);
		Q=.5*RHO*VT^2;
		GAMT=atan2(-VT2,VT1);
		AT1=-32.2*Q*cos(GAMT)/BETA;
		AT2=-32.2+32.2*Q*sin(GAMT)/BETA;
		FLAG=1;
	end;
	FLAG=0;
	RT1=.5*(RT1OLD+RT1+H*VT1);
	RT2=.5*(RT2OLD+RT2+H*VT2);
	VT1=.5*(VT1OLD+VT1+H*AT1);
	VT2=.5*(VT2OLD+VT2+H*AT2);
	S=S+H;
 	if S>=.09999
		S=0.;
		ATG=sqrt(AT1^2+AT2^2)/32.2;
		RT1K=RT1/1000.;
		RT2K=RT2/1000.;
		VT=sqrt(VT1^2+VT2^2);
		count=count+1;
		ArrayT(count)=T;
		ArrayRT1K(count)=RT1K;
		ArrayRT2K(count)=RT2K;
		ArrayATG(count)=ATG;
		ArrayVT(count)=VT;
	end
end
figure
plot(ArrayRT1K,ArrayRT2K),grid
xlabel('Downrange (Kft)')
ylabel('Altitude (Kft) ')
figure
plot(ArrayRT2K,ArrayATG),grid
xlabel('Altitude (Kft)')
ylabel('Acceleration (G) ')
figure
plot(ArrayRT2K,ArrayVT),grid
xlabel('Altitude (Kft)')
ylabel('Velocity (Ft/Sec) ')
clc
output=[ArrayT',ArrayRT1K',ArrayRT2K',ArrayATG',ArrayVT'];
save datfil.txt output -ascii
disp 'simulation finished'
