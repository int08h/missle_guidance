% Preallocation
clear
H=zeros(1,10000);
X=zeros(1,10000);
count=0;
XMAX=6;
XMIN=-6;
RANGE=XMAX-XMIN;
TMP=1./sqrt(6.28);
N=100;
BIN=50;
for I=1:N
	SUM=0;
	for J=1:12
		RAND=rand(1);
		SUM=SUM+RAND;
	end;
   	X(I)=SUM-6;
end;
for I=1:BIN
	H(I)=0;
end;
% FIX	Round towards zero.
% 	FIX(X) rounds the elements of X to the nearest integers
%  	towards zero.
for I=1:N
        K=fix(((X(I)-XMIN)/RANGE)*BIN)+.99;
        if K < 1, K=1; end;
        if K > BIN, K=BIN; end;
%CORRECTION HERE >>>
	K=round(K);
        H(K)=H(K)+1;
end;
for K=1:BIN
	PDF=(H(K)/N)*BIN/RANGE;
   	AB=XMIN+K*RANGE/BIN;
   	TH=TMP*exp(-AB*AB/2.);
	count=count+1;
	ArrayAB(count)=AB;
	ArrayPDF(count)=PDF;
	ArrayTH(count)=TH;
end;
figure
plot(ArrayAB,ArrayPDF,ArrayAB,ArrayTH),grid
title('Sample Gaussian distribution')
xlabel('X')
ylabel('Probability Density Function')
clc
output=[ArrayAB',ArrayPDF',ArrayTH'];
save datfil output  -ascii
disp 'simulation finished'

