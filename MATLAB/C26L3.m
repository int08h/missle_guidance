clear
count=0;
ORDER=4;
WACT=100.;
ZACT=.65;
KDEL=9000.;;
WR=2.;
DELCMX=30.;
PHIMX=10.;
PHIDMX=300.;
G=zeros([ORDER,1]);
A=zeros([ORDER,ORDER]);
F=zeros([ORDER,ORDER]);
S=zeros([ORDER,ORDER]);
G(4,1)=WACT*WACT;
F(1,2)=1;
F(2,2)=-WR;
F(2,3)=KDEL;
F(3,4)=1.;
F(4,3)=-WACT*WACT;
F(4,4)=-2.*ZACT*WACT;
A(1,1)=(DELCMX/PHIMX)^2;
A(2,2)=(DELCMX/PHIDMX)^2;
T=0;
H=.0002;
S1=0;
TF=5.;
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
			SD=S*F+F'*S-S*G*G'*S+A;
			C=GT*S;
			FLAG=1;
		end
		FLAG=0;
		H2=.5*H;
		HSDP=H2*SD;
		SS=SOLD+S;
		SSP=.5*SS;
		S=SSP+HSDP;
		if S1>=.0499999
			S1=0;
			C1=C(1,1);
			C2=C(1,2);
			C3=C(1,3);
			C4=C(1,4);
			count=count+1;
			ArrayT(count)=T;
			ArrayC1(count)=C1;
			ArrayC2(count)=C2;
			ArrayC3(count)=C3;
			ArrayC4(count)=C4;
		 end
end
output=[ArrayT',ArrayC1',ArrayC2',ArrayC3',ArrayC4'];
save datfil.txt output -ascii
disp 'simulation finished'
clc
figure
plot(ArrayT,ArrayC1),grid
xlabel('Time (s) ')
ylabel('C1')
figure
plot(ArrayT,ArrayC2),grid
xlabel('Time (s) ')
ylabel('C2')
figure
plot(ArrayT,ArrayC3),grid
xlabel('Time (s) ')
ylabel('C3')
figure
plot(ArrayT,ArrayC4),grid
xlabel('Time (s) ')
ylabel('C4')
C1
C2
C3
C4
		
		

				 
