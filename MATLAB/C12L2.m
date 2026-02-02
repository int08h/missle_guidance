clear all
close all
ITERM=1;
G=32.2;
SIGNOISE=25.;
X=200000.;
XD=-6000.;
BETA=500.;
XH=200025.;
XDH=-6150.;
BETAH=800.;
BETAINV=1./BETA;
BETAINVH=1./BETAH;
ORDER=3;
TS=.1;
TF=30.;
PHIS=0.;
T=0.;
S=0.;
H=.001;
HP=.001;
PHI=zeros(ORDER,ORDER);
P=[SIGNOISE*SIGNOISE 0 0
   ;0 20000. 0;
   0 0 (BETAINV-BETAINVH)^2];
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
		F21=-G*RHOH*XDH*XDH*BETAINVH/44000.;
		F22=RHOH*G*XDH*BETAINVH;
		F23=.5*RHOH*XDH*XDH*G;
		F=[0  1  0;
		   F21  F22  F23;
		   0  0  0];
		PHI=IDNP+F*TS;
		Q=[0  0  0;
		   0  F23*F23*PHIS*TS*TS*TS/3  F23*PHIS*TS*TS/2;
		   0  F23*PHIS*TS*TS/2  PHIS*TS];
		M=PHI*P*PHI'+Q;
		K = M*HMAT'/(HMAT*M*HMAT' + RMAT);
		P = (IDNP - K*HMAT)*M;
		XNOISE=SIGNOISE*randn;
		BETAH=1./BETAINVH;
		[XB,XDB]=PROJECTC12L2(T,TS,XH,XDH,BETAINVH,HP);
		RES=X+XNOISE-XB;
		XH=XB+K(1,1)*RES;
		XDH=XDB+K(2,1)*RES;
		BETAINVH=BETAINVH+K(3,1)*RES;
		ERRX=X-XH;
		SP11=sqrt(P(1,1));
		ERRXD=XD-XDH;
		SP22=sqrt(P(2,2));
		ERRBETAINV=1./BETA-BETAINVH;
		SP33=sqrt(P(3,3));
		BETAH=1./BETAINVH;
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
		ArrayERRBETAINV(count)=ERRBETAINV;
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
plot(ArrayT,ArrayERRBETAINV,ArrayT,ArraySP33,ArrayT,ArraySP33P),grid
xlabel('Time (Sec)')
ylabel('Error in Estimate of 1/BETA (Ft^2/Lb)')
axis([0 30 -.0008 .0008])
clc
output=[ArrayT',ArrayX',ArrayXH',ArrayXD',ArrayXDH',ArrayBETA',ArrayBETAH'];
save datfil.txt output  -ascii
output=[ArrayT',ArrayERRX',ArraySP11',ArraySP11P',ArrayERRXD',ArraySP22',...
ArraySP22P',ArrayERRBETAINV',ArraySP33',ArraySP33P'];
save covfil.txt output  -ascii
disp 'simulation finished'

