clear
n=0;		
VC=5.*3280.;
XNTIC=161.;
VM=3000.;
HEDEG=20.;
SIGNOISE=10.;
TS=.01;
for TF=.2:.2:10.0,
	TS2=TS*TS;
	TS3=TS2*TS;
	TS4=TS3*TS;
	TS5=TS4*TS;
	PHIS=XNTIC*XNTIC/TF;
	PHIN=SIGNOISE*SIGNOISE*TS;
	SIGPOS=SIGNOISE;
	SIGN2=SIGPOS^2;
	P11=SIGN2;
	P12=0.;
	P13=0.;
	P22=(VM*HEDEG/57.3)^2;
	P23=0.;
	P33=XNTIC*XNTIC;
	T=0.;
	for T=TS:TS:TF,
		TGO=TF-T+.000001;
		RTM=VC*TGO;
		SIGPOS=SIGNOISE;
		SIGN2=SIGPOS^2;
		M11=P11+TS*P12+.5*TS2*P13+TS*(P12+TS*P22+.5*TS2*P23);
		M11=M11+.5*TS2*(P13+TS*P23+.5*TS2*P33)+TS5*PHIS/20.;
		M12=P12+TS*P22+.5*TS2*P23+TS*(P13+TS*P23+.5*TS2*P33);
     		M12=M12+TS4*PHIS/8.;
		M13=P13+TS*P23+.5*TS2*P33+PHIS*TS3/6.;
		M22=P22+TS*P23+TS*(P23+TS*P33)+PHIS*TS3/3.;
		M23=P23+TS*P33+.5*TS2*PHIS;
		M33=P33+PHIS*TS;
		K1=M11/(M11+SIGN2);
		K2=M12/(M11+SIGN2);
		K3=M13/(M11+SIGN2);
		P11=(1.-K1)*M11;
		P12=(1.-K1)*M12;
		P13=(1.-K1)*M13;
		P22=-K2*M12+M22;
		P23=-K2*M13+M23;
		P33=-K3*M13+M33;
		SP11=sqrt(P11);
	end;
 	FORM=sqrt(2.*(PHIS^.16667)*(PHIN^.83333));
 	n=n+1;
	ArrayT(n)=T;
	ArrayFORM(n)=FORM;
	ArraySP11(n)=SP11;
end;
figure
plot(ArrayT',ArrayFORM',ArrayT',ArraySP11'),grid
title('RMS miss for various flight times')
xlabel('Flight Time (S)')
ylabel('RMS MISS (Ft) ')
clc
output=[ArrayT',ArrayFORM',ArraySP11'];
save datfil.txt output -ascii
disp('Simulation Complete')
