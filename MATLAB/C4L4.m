clear
TAU=0.2;
PHI=1.0;
count=0;
T=0;
H=.01;
SIG=sqrt(PHI/H);
Y=0;
while T <= 5.0
	X=SIG*randn;
   	YOLD=Y;
   	STEP=1;
	FLAG=0;
	while STEP <=1
		if FLAG ==1
			STEP=2;
			Y=Y+H*YD;
         		T=T+H;
		end
		YD=(X-Y)/TAU;
		FLAG=1;
	end
	FLAG=0;
   	Y=(YOLD+Y)/2.+.5*H*YD;
   	SIGPLUS=sqrt(PHI*(1.-exp(-2.*T/TAU))/(2.*TAU));
   	SIGMINUS=-SIGPLUS;
	count=count+1;
	ArrayT(count)=T;
	ArrayY(count)=Y;
	ArraySIGPLUS(count)=SIGPLUS;
	ArraySIGMINUS(count)=SIGMINUS;
end
figure
plot(ArrayT,ArrayY,ArrayT,ArraySIGPLUS,ArrayT,ArraySIGMINUS),grid
title('Simulation of low-pass filter driven by white noise')
xlabel('Time (S)')
ylabel('Y')
clc
output=[ArrayT',ArrayY',ArraySIGPLUS',ArraySIGMINUS'];
save datfil.txt output  -ascii
disp 'simulation finished'
