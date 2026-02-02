clear
count=0;
for I=2:160
	W=10^(.025*I-1);
	TOP1=13500*sqrt(1+(W/29.1)^2);
	TOP2=sqrt((1-W*W/(173.*173.))^2+(2.*.74*W/173)^2);
	BOT1=W*sqrt(1+(W/2.)^2);
	BOT2=sqrt((1-W*W/(100.*100.))^2+(2.*.65*W/100)^2);
	XMAG=TOP1*TOP2/(BOT1*BOT2);
	GAIN=20*log10(XMAG);
	count=count+1;
	ArrayW(count)=W;
	ArrayGAIN(count)=GAIN;
end
figure
semilogx(ArrayW,ArrayGAIN),grid
xlabel('Frequency (Rad/Sec)')
ylabel('Gain (Db)')
clc
output=[ArrayW',ArrayGAIN'];
save datfil.txt output -ascii
 	