clear
W=1;
GAM=.00001;
TF=10;
F=zeros([4,4]);
S=zeros([4,4]);
count=0;
G(1,1)=0;
G(2,1)=-1.;
G(3,1)=0;
G(4,1)=0.;
F(1,2)=1;
F(2,3)=1;
F(3,4)=1.;
F(4,3)=-W*W;
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
		NP=C2*T;
		XNPP=3.;
		C1TH=XNPP/(T*T);
		C2TH=XNPP/T;
		C3TH=XNPP*((1.-cos(W*T))/W^2)/T^2;
		C4TH=XNPP*((W*T-sin(W*T))/W^3)/T^2;
		count=count+1;
		ArrayT(count)=T;
		ArrayC1(count)=C1;
		ArrayC1TH(count)=C1TH;
		ArrayC2(count)=C2;
		ArrayC2TH(count)=C2TH;
		ArrayC3(count)=C3;
		ArrayC3TH(count)=C3TH;
		ArrayC4(count)=C4;
		ArrayC4TH(count)=C4TH;
		ArrayNP(count)=NP;
		ArrayXNPP(count)=XNPP;
	end
end
output=[ArrayT',ArrayC1',ArrayC1TH',ArrayC2',ArrayC2TH',...
	ArrayC3',ArrayC3TH',ArrayC4',ArrayC4TH',ArrayNP',ArrayXNPP'];
save datfil.txt output -ascii
disp 'simulation finished'
clc
figure
plot(ArrayT,ArrayC3,ArrayT,ArrayC3TH),grid
xlabel('Time (s) ')
ylabel('C3')	
axis([0 10 0 1.5])
figure
plot(ArrayT,ArrayC4,ArrayT,ArrayC4TH),grid
xlabel('Time (s) ')
ylabel('C4')	
axis([0 10 0 1])



