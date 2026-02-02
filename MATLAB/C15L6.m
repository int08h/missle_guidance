clear
n=0;
GM=1.4077E16;
A=2.0926E7;
CONST=sqrt(GM/A);
for DISTKM=100:100:20000,
	PHI=DISTKM*3280./A ;
	GAM=3.14159/4.-PHI/4.;
	GAMDEG=57.3*GAM;
	TOP=GM*(1.-cos(PHI)); 
	TEMP=A*cos(GAM)/A-cos(PHI+GAM); 
	BOT=A*cos(GAM)*TEMP; 
	V=sqrt(TOP/BOT); 
	XLAM=A*V*V/GM;
	TOP1=tan(GAM)*(1-cos(PHI))+(1-XLAM)*sin(PHI); 
	BOT1P=(1-cos(PHI))/(XLAM*cos(GAM)*cos(GAM)); 
	BOT1=(2-XLAM)*(BOT1P+cos(GAM+PHI)/cos(GAM)); 
	TOP2=2*cos(GAM); 
	BOT2=XLAM*((2/XLAM-1)^1.5); 
	TOP3=sqrt(2/XLAM-1); 
	BOT3=cos(GAM)/tan(PHI/2)-sin(GAM); 
	TEMP=(TOP2/BOT2)*atan2(TOP3,BOT3); 
	TF=A*(TOP1/BOT1+TEMP)/(V*cos(GAM));
	n=n+1;
	ArrayDISTKM(n)=DISTKM;
	ArrayTF(n)=TF;
end
figure
plot(ArrayDISTKM,ArrayTF),grid
xlabel('Distance (km)')
ylabel('Flight Time (s)')
clc
output=[ArrayDISTKM',ArrayTF'];
save datfil.txt output -ascii
disp 'simulation finished'

