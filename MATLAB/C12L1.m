clear all
close all
ITERM=1;
SIGNOISE=25.;
X=200000.;
XD=-6000.;
BETA=500.;
XH=200025.;
XDH=-6150.;
BETAH=800.;
ORDER=3;
TS=.1;
TF=30.;
PHIS=0.;
T=0.;
S=0.;
H=.001;
HP=.001;
P=[SIGNOISE*SIGNOISE 0 0;
   0 20000. 0;
   0 0 300.^2];
IDNP=eye(ORDER); 	
HMAT=[1 0 0];
RMAT=SIGNOISE^2;
count=0;
while T<=TF 
	XOLD=X;
	XDOLD=XD;
 	XDD=.0034*32.2*XD*XD*exp(-X/22000.)/(2.*BETA)-32.2;
 	X=X+H*XD;
	XD=XD+H*XDD;
 	T=T+H;
	XDD=.0034*32.2*XD*XD*exp(-X/22000.)/(2.*BETA)-32.2;
 	X=.5*(XOLD+X+H*XD);
	XD=.5*(XDOLD+XD+H*XDD);
 	S=S+H;
	if S>=(TS-.00001)
		S=0.;
		RHOH=.0034*exp(-XH/22000.);
		F21=-32.2*RHOH*XDH*XDH/(44000.*BETAH);
		F22=RHOH*32.2*XDH/BETAH;
		F23=-RHOH*32.2*XDH*XDH/(2.*BETAH*BETAH);
		F=[0  1  0;
		   F21  F22  F23;
		   0  0  0 ];
		if ITERM==1
			PHI=IDNP+F*TS;
		else
			PHI=IDNP+F*TS+F*F*TS*tTS/2;
		end
		Q=[0  0  0;
			0  F23*F23*PHIS*TS*TS*TS/3  F23*PHIS*TS*TS/2;
			0  F23*PHIS*TS*TS/2  PHIS*TS];
		
		M=PHI*P*PHI'+Q;
		K = M*HMAT'/(HMAT*M*HMAT' + RMAT);
		P = (IDNP - K*HMAT)*M;
 		XNOISE=SIGNOISE*randn;
 		[XB,XDB]=PROJECTC12L1(T,TS,XH,XDH,BETAH,HP);
		RES=X+XNOISE-XB;
		XH=XB+K(1,1)*RES;
		XDH=XDB+K(2,1)*RES;
		BETAH=BETAH+K(3,1)*RES;
		ERRX=X-XH;
		SP11=sqrt(P(1,1));
		ERRXD=XD-XDH;
		SP22=sqrt(P(2,2));
		ERRBETA=BETA-BETAH;
		SP33=sqrt(P(3,3));
		SP11P=-SP11;
		SP22P=-SP22;
		SP33P=-SP33;
		count=count+1;
		ArrayT(count)=T;
		ArrayX(count)=X;
		ArrayXH(count)=XH;
		ArrayXD(count)=XD;
		ArrayXDH(count)=XDH;
		ArrayBETA(count)=BETA;
		ArrayBETAH(count)=BETAH;
		ArrayERRX(count)=ERRX;
		ArraySP11(count)=SP11;
		ArraySP11P(count)=SP11P;
		ArrayERRXD(count)=ERRXD;
		ArraySP22(count)=SP22;
		ArraySP22P(count)=SP22P;
		ArrayERRBETA(count)=ERRBETA;
		ArraySP33(count)=SP33;
		ArraySP33P(count)=SP33P;
	end
end
figure
plot(ArrayT,ArrayERRX,ArrayT,ArraySP11,ArrayT,ArraySP11P),grid
xlabel('Time (Sec)')
ylabel('Error in Estimate of Altitude (Ft)')
axis([0 30 -25 25])
figure
plot(ArrayT,ArrayERRXD,ArrayT,ArraySP22,ArrayT,ArraySP22P),grid
xlabel('Time (Sec)')
ylabel('Error in Estimate of Velocity (Ft/Sec)')
axis([0 30 -25 25])
figure
plot(ArrayT,ArrayERRBETA,ArrayT,ArraySP33,ArrayT,ArraySP33P),grid
xlabel('Time (Sec)')
ylabel('Error in Estimate of Ballistic Coefficient (Lb/Ft^2)')
clc
output=[ArrayT',ArrayX',ArrayXH',ArrayXD',ArrayXDH',ArrayBETA',ArrayBETAH'];
save datfil.txt output  -ascii
output=[ArrayT',ArrayERRX',ArraySP11',ArraySP11P',ArrayERRXD',ArraySP22',...
ArraySP22P',ArrayERRBETA',ArraySP33',ArraySP33P'];
save covfil.txt output  -ascii
disp 'simulation finished'
		
