clear
count=0;
TAU=1.;
GAM=.00001;
T=0.;
H=.01;
S=0.;
X1=0.;
X2=0.;
X3=0.;
while ~(T >= (10.-.0001))
	S=S+H;
	X1OLD=X1;
	X2OLD=X2;
	X3OLD=X3;
	STEP=1;
	FLAG=0;
	while STEP<=1
		if FLAG==1 
         		STEP=2;
			X1=X1+H*X1D;
			X2=X2+H*X2D;
			X3=X3+H*X3D;
			T=T+H;
		end
		X1D=(T-X1)/TAU;
		X2D=X1^2;
		D=X2+GAM;
		X3D=-X3/TAU+T;
		FLAG=1;
	end
	FLAG=0;
  	X1=.5*(X1OLD+X1+H*X1D);
	X2=.5*(X2OLD+X2+H*X2D);
	X3=.5*(X3OLD+X3+H*X3D);
	if S>=.09999
		S=0.;
		PZ=X1/D;
		XNP=PZ*T*T;
		C1=PZ;
		C2=PZ*T;
		C3=.5*PZ*T*T;
		C4=-X3*PZ;
		XS=T/TAU;
		TOP=6.*XS*XS*(exp(-XS)-1.+XS);
		BOT1=2*XS*XS*XS+3.+6.*XS-6.*XS*XS;
		BOT2=-12.*XS*exp(-XS)-3.*exp(-2.*XS);
		XNPTH=TOP/(BOT1+BOT2+.0001);
		C4TH=-XNPTH*(exp(-XS)+XS-1.)/(XS*XS);
		count=count+1;
		ArrayT(count)=T;
		ArrayC4(count)=C4;
		ArrayC4TH(count)=C4TH;
		ArrayXNP(count)=XNP;
		ArrayXNPTH(count)=XNPTH;
	end
end
output=[ArrayT',ArrayC4',ArrayC4TH',ArrayXNP',ArrayXNPTH'];
save datfil.txt output -ascii
disp 'simulation finished'
clc
figure
semilogy(ArrayT,ArrayXNP,ArrayT,ArrayXNPTH),grid
xlabel('Time (s) ')
ylabel('NP')	
axis([0 10 .1 40])
figure
semilogy(ArrayT,-ArrayC4,ArrayT,-ArrayC4TH),grid
xlabel('Time (s) ')
ylabel('-C4')	
axis([0 10 .1 20])

