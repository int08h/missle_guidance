clear
F=zeros([4,4]);
X=zeros([4,4]);
XOLD=zeros([4,4]);
Q=zeros([4,4]);
K0=zeros([4,4]);
K1=zeros([4,4]);
K2=zeros([4,4]);
K3=zeros([4,4]);
XD=zeros([4,4]);
A=zeros([1,4]);
AXAT=zeros([1,1]);
count=0;
T=0.;
TNEW=T;
S=0.;
H=.01;
XNP=3.;
TAU=1.;
XNT=96.6;
VC=4000.;
TF=10.;
TGO=TF-T+.00001;
PHIS=XNT*XNT/TF;
F(1,2)=1.;
F(2,1)=-XNP/(TAU*TGO);
F(2,3)=1.;
F(2,4)=XNP*VC/TAU;
F(4,1)=1./(TAU*VC*TGO);
F(4,4)=-1./TAU;
Q(3,3)=PHIS;
while ~(T >= (TF-.0001))
	S=S+H;
	XOLD=X;
	STEP=1;
	while ( (STEP == 1) | (STEP == 2) | (STEP == 3))
		TGO=TF-TNEW+.00001;
		F(2,1)=-XNP/(TAU*TGO);
 		F(4,1)=1./(TAU*VC*TGO);
		XD=(F*X)+(F*X)'+Q;
		if (STEP == 1)
			STEP=2;
			K0=XD;
 			TNEW=T+.5*H;
			X=XOLD+.5.*H.*K0;
		elseif (STEP ==2)
			STEP=3;
			K1=XD;
 			TNEW=T+.5*H;
			X=XOLD+.5.*H.*K1;
		else
			STEP=4;
			K2=XD;
		 	TNEW=T+H;
			X=XOLD+H.*K2;
		end
	end
	TGO=TF-TNEW+.00001;
	F(2,1)=-XNP/(TAU*TGO);
 	F(4,1)=1./(TAU*VC*TGO);
 	XD=(F*X)+(F*X)'+Q;
 	K3=XD;
 	T=TNEW;
	X=XOLD+H.*(K0+2.*(K1+K2)+K3)./6;
	if S>=.09999
		S=0.;
		A(1,1)=XNP/(TAU*TGO);
		A(1,2)=0.;
		A(1,3)=0.;
		A(1,4)=-XNP*VC/TAU;
		AXAT=A*X*A';
		SIGY=sqrt(X(1,1));
		SIGNL=sqrt(AXAT(1,1));
		count=count+1;
		ArrayT(count)=T;
		ArraySIGY(count)=SIGY;
		ArraySIGNLG(count)=SIGNL/32.2;
	end
end
figure
plot(ArrayT,ArraySIGY)
xlabel('Time (s) ')
ylabel('Standard Deviation of Relative Position (Ft)')
figure
plot(ArrayT,ArraySIGNLG)
xlabel('Time (s) ')
ylabel('Standard Deviation of Acceleration (G)')
axis([0 10 0 20])
clc
output=[ArrayT', ArraySIGY',ArraySIGNLG'];
save datfil.txt output -ascii
disp 'simulation finished'
SIGY

