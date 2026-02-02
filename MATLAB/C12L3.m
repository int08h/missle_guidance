clear all
close all
ITERM=1;
G=32.2;
SIGNOISE=25.;
RMAT=SIGNOISE^2;
X=200000.;
XD=-6000.;
BETA=500.;
XH=200025.;
XDH=-6150.;
XDDH=0.;
XNT=322.;
ORDER=3;
TS=.1;
TF=30;
PHIS=XNT*XNT/TF;
T=0.;
S=0.;
H=.001;
HP=.001;
TS2=TS*TS;
TS3=TS2*TS;
TS4=TS3*TS;
TS5=TS4*TS;
PHI=[1 TS .5*TS*TS;
	 0 1 TS; 
	 0 0 1];
P=[SIGNOISE*SIGNOISE 0 0;
   0 20000 0;
   0 0 XNT^2];
IDNP=eye(ORDER);
Q=[TS5*PHIS/20  TS4*PHIS/8  PHIS*TS3/6;
   TS4*PHIS/8  PHIS*TS3/3  .5*TS2*PHIS;
   PHIS*TS3/6  .5*TS2*PHIS  PHIS*TS];
HMAT=[1 0 0];
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
		M=PHI*P*PHI'+Q;
		K = M*HMAT'/(HMAT*M*HMAT' + RMAT);
		P = (IDNP - K*HMAT)*M;		XNOISE=SIGNOISE*randn;
		RES=X+XNOISE-XH-TS*XDH-.5*TS*TS*XDDH;
		XH=XH+TS*XDH+.5*TS*TS*XDDH+K(1,1)*RES;
		XDH=XDH+TS*XDDH+K(2,1)*RES;
		XDDH=XDDH+K(3,1)*RES;
		RHOH=.0034*exp(-XH/22000.);
		BETAH=16.1*RHOH*XDH*XDH/(XDDH+32.2);
		ERRX=X-XH;
		SP11=sqrt(P(1,1));
		ERRXD=XD-XDH;
		SP22=sqrt(P(2,2));
		ERRXDD=XDD-XDDH;
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
		ArrayXDD(count)=XDD;
		ArrayXDDH(count)=XDDH;
		ArrayERRX(count)=ERRX;
		ArraySP11(count)=SP11;
		ArraySP11P(count)=SP11P;
		ArrayERRXD(count)=ERRXD;
		ArraySP22(count)=SP22;
		ArraySP22P(count)=SP22P;
		ArrayERRXDD(count)=ERRXDD;
		ArraySP33(count)=SP33;
		ArraySP33P(count)=SP33P;
		end
end
figure
plot(ArrayT,ArrayERRX,ArrayT,ArraySP11,ArrayT,ArraySP11P),grid
xlabel('Time (Sec)')
ylabel('Error in Estimate of Altitude (Ft)')
axis([0 30 -50 50])
figure
plot(ArrayT,ArrayERRXD,ArrayT,ArraySP22,ArrayT,ArraySP22P),grid
xlabel('Time (Sec)')
ylabel('Error in Estimate of Velocity (Ft/Sec)')
axis([0 30 -150 150])
figure
plot(ArrayT,ArrayERRXDD,ArrayT,ArraySP33,ArrayT,ArraySP33P),grid
xlabel('Time (Sec)')
ylabel('Error in Estimate of Acceleration (Ft/Sec^2)')
axis([0 30 -200 200])
clc
output=[ArrayT',ArrayX',ArrayXH',ArrayXD',ArrayXDH',ArrayXDD',ArrayXDDH'];
save datfil.txt output  -ascii
output=[ArrayT',ArrayERRX',ArraySP11',ArraySP11P',ArrayERRXD',ArraySP22',...
ArraySP22P',ArrayERRXDD',ArraySP33',ArraySP33P'];
save covfil.txt output  -ascii
disp 'simulation finished'

	
	
