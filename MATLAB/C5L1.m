% Preallocation
clear
K0=zeros([1,6]);
K1=zeros([1,6]);
K2=zeros([1,6]);
K3=zeros([1,6]);
count=0;
XIN=1.;
W=20.;
Y=0.;
YD=0.;
T=0.;
H=.01;
S=0.;
while ~(T >= 1.)
	S=S+H;
	YOLD=Y;
	YDOLD=YD;
	STEP=1;
	while ( (STEP == 1) | (STEP == 2) | (STEP == 3))
		YDD=W*XIN-W*W*Y;
		if (STEP == 1)
			STEP=2;
			K0(1,1)=YD;
			K0(1,2)=YDD;
 			TNEW=T+.5*H;
			Y=YOLD+.5.*H.*K0(1,1);
			YD=YDOLD+.5.*H.*K0(1,2);
		elseif (STEP ==2)
			STEP=3;
			K1(1,1)=YD;
			K1(1,2)=YDD;
 			TNEW=T+.5*H;
			Y=YOLD+.5.*H.*K1(1,1);
			YD=YDOLD+.5.*H.*K1(1,2);
		else
			STEP=4;
			K2(1,1)=YD;
			K2(1,2)=YDD;
		 	TNEW=T+H;
			Y=YOLD+H.*K2(1,1);
			YD=YDOLD+H.*K2(1,2);
		end
	end
	YDD=W*XIN-W*W*Y;
	K3(1,1)=YD;
	K3(1,2)=YDD;
 	T=TNEW;
	Y=YOLD+H.*(K0(1,1)+2.*(K1(1,1)+K2(1,1))+K3(1,1))./6;
	YD=YDOLD+H.*(K0(1,2)+2.*(K1(1,2)+K2(1,2))+K3(1,2))./6;
	count=count+1;
	ArrayT(count)=T;
  	ArrayY(count)=Y;
end
figure
plot(ArrayT,ArrayY)
xlabel('Time (s) ')
ylabel('Y')
title('Forth-order Runge-Kutta: Second Order Network')
output=[ArrayT', ArrayY'];
save datfil.txt output -ascii
disp 'simulation finished'
