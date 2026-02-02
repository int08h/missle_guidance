clear
count=0;
N=100;
for I=1:N       
  	SUM=0;
   	for J=1:12
      		RAND=rand(1);
		SUM=SUM+RAND;
	end;
   	X=SUM-6;
	count=count+1;
	ArrayI(count)=I;
	ArrayX(count)=X;
end;
figure
plot(ArrayI,ArrayX),grid
title('One hundred random numbers with Gaussian distribution')
xlabel('Number')
ylabel('Value')
clc
output=[ArrayI',ArrayX'];
save datfil output  -ascii
disp 'simulation finished'

