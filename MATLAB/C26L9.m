clear
count=0;
for I=2:160
	W=10^(.025*I-1);
	TOP1=13500*sqrt(1+(W/29.)^2);
	TOP2=sqrt((1-W*W/(160.*160.))^2+(2.*.83*W/160)^2);
	TOP3=sqrt((1-W*W/(216.*216.))^2+(2.*.45*W/216)^2);
	BOT1=W*sqrt(1+(W/2.)^2);
	BOT2=sqrt((1-W*W/(100.*100.))^2+(2.*.65*W/100)^2);
	BOT3=sqrt((1-W*W/(200.*200.))^2+(2.*.5*W/200)^2);
	XMAG=TOP1*TOP2*TOP3/(BOT1*BOT2*BOT3);
	GAIN=20*log10(XMAG);
	count=count+1;
	ArrayW(count)=W;
	ArrayGAIN(count)=GAIN;
end
figure
semilogx(ArrayW,ArrayGAIN),grid
xlabel('Frequency (Rad/Sec)')
ylabel('Gain (Db)')
axis([.1 1000 -60 40])
figure
semilogx(ArrayW,ArrayPHASE),grid
xlabel('Frequency (Rad/Sec)')
ylabel('Phase (Deg)')
axis([.1 1000 -400 100])
clc
output=[ArrayW',ArrayGAIN',ArrayPHASE'];
save datfil.txt output -ascii
disp 'simulation finished' 	
 	