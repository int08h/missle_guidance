clear all
close all
count=0;
W=2;
T=0;
S=0;
X=0;
XD=W;
H=.01;
while T<=10
	S=S+H;	
	XOLD=X;
	XDOLD=XD;
	XDD=-W*W*X;
	X=X+H*XD;
	XD=XD+H*XDD;
	T=T+H;
	XDD=-W*W*X;
	X=.5*(XOLD+X+H*XD);
	XD=.5*(XDOLD+XD+H*XDD);
	if S>=.09999
		S=0;
		XTHEORY=sin(W*T);
		count=count+1;
      		ArrayT(count)=T;
      		ArrayX(count)=X;
      		ArrayXTHEORY(count)=XTHEORY;
	end
end
figure
plot(ArrayT,ArrayX,ArrayT,ArrayXTHEORY),grid
xlabel('Time (Sec)')
ylabel('X&XTHEORY ')
axis([0 10 -2 2])
clc
output=[ArrayT',ArrayX',ArrayXTHEORY'];
save datfil.txt output -ascii
disp 'simulation finished'
