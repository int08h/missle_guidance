clear
TAU=1.;
WZ=5.;
W=20.;
Z=.7;
TF=10.;
GAM=.00001;
F=zeros([6,6]);
S=zeros([6,6]);
count=0;
G(1,1)=0.;
G(2,1)=0.;
G(3,1)=0.;
G(4,1)=0.;
G(5,1)=0.;
G(6,1)=W*W/TAU;
F(1,2)=1;
F(2,3)=1;
F(2,4)=-1.;
F(2,6)=1./WZ^2;
F(4,5)=1.;
F(5,6)=1.;
F(6,4)=-W*W/TAU;
F(6,5)=-W*W*(2.*Z/W+TAU)/TAU;
F(6,6)=-W*W*(1./W^2+2.*Z*TAU/W)/TAU;
S(1,1)=1;
T=0;
H=.0001;
S1=0;
while ~(T >= (TF-.0001))
	S1=S1+H;
	SOLD=S;
	STEP=1;
 	FLAG=0;
	while STEP<=1
		if FLAG==1 
         		STEP=2;
         		HSD=H*SD;
         		S=S+HSD;
 			T=T+H;
 		end
 		SF=S*F;
 		GT=G';
 		GTS=GT*S;
 		GAMINV=1./GAM;
 		C=GAMINV*GTS;
 		SFT=SF';
 		CT=C';
 		CTC=CT*C;
 		CTBC=GAM*CTC;
 		SFSFT=SF+SFT;
		SD=SFSFT-CTBC;
		FLAG=1;
	end
	FLAG=0;
	H2=.5*H;
	HSDP=H2*SD;
	SS=SOLD+S;
	SSP=.5*SS;
	S=SSP+HSDP;
 	if S1>=.009999
		S1=0;
		C1=-C(1,1);
		C2=-C(1,2);
		C3=-C(1,3);
		C4=-C(1,4);
		C5=-C(1,5);
		C6=-C(1,6);
		NP=C2*T;
		count=count+1;
		ArrayT(count)=T;
		ArrayC4(count)=C4;
		ArrayC5(count)=C5;
		ArrayC6(count)=C6;
		ArrayNP(count)=NP;
	end
end
output=[ArrayT',ArrayC4',ArrayC5',ArrayC6',ArrayNP'];
save datfil.txt output -ascii
disp 'simulation finished'
clc
figure
plot(ArrayT,ArrayNP),grid
xlabel('Time (s) ')
ylabel('NP')	
axis([0 10 -10 50])
