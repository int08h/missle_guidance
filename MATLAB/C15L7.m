clear
n=0;
GM=1.4077E16;
A=2.0926E7;
CONST=sqrt(GM/A);
DISTKM=10000.;
PHI=DISTKM*3280./A;
%for GAMDEG=0:1:70
for GAMDEG=160:1:180
	GAM=GAMDEG/57.3;
	TOP=GM*(1.-cos(PHI));
	TEMP=A*cos(GAM)/A-cos(PHI+GAM);
	BOT=A*cos(GAM)*TEMP;
	V=sqrt(TOP/BOT);
	VKM=V/3280.;
	XLAM=A*V*V/GM;
	TOP1=tan(GAM)*(1-cos(PHI))+(1-XLAM)*sin(PHI) ;
	BOT1P=(1-cos(PHI))/(XLAM*cos(GAM)*cos(GAM));;
	BOT1=(2-XLAM)*(BOT1P+cos(GAM+PHI)/cos(GAM));
	TOP2=2*cos(GAM); 
	BOT2=XLAM*((2/XLAM-1)^1.5) ;
	TOP3=sqrt(2/XLAM-1) ;
	BOT3=cos(GAM)/tan(PHI/2)-sin(GAM); 
	TEMP=(TOP2/BOT2)*atan2(TOP3,BOT3); 
	TF=A*(TOP1/BOT1+TEMP)/(V*cos(GAM)); 
	n=n+1;
	ArrayGAMDEG(n)=GAMDEG;
	ArrayVKM(n)=VKM;
	ArrayTF(n)=TF;
	ArrayXLAM(n)=XLAM;
end
figure
plot(ArrayGAMDEG,ArrayVKM),grid
xlabel('Flight Path Angle (deg)')
ylabel('Velocity (km/s)')
figure
plot(ArrayGAMDEG,ArrayTF),grid
xlabel('Flight Path Angle (deg)')
ylabel('Flight Time (s)')
clc
output=[ArrayGAMDEG',ArrayVKM',ArrayTF',ArrayXLAM'];
save datfil.txt output -ascii
disp 'simulation finished'

