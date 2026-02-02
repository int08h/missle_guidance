clear
G=.5;
X=1.;
TS=.1;
Y=0.;
T=0.;
N=0;
count=0;
YTHEORY=1.-(1.-G)^N;
for N=1:20
   Y=Y+G*(X-Y);
   T=N*TS;
   YTHEORY=1.-(1.-G)^N;
   count=count+1;
   ArrayT(count)=T;
   ArrayY(count)=Y;
   ArrayYTHEORY(count)=YTHEORY;
end;
figure
plot(ArrayT,ArrayY,ArrayT,ArrayYTHEORY),grid
title('Output')
xlabel('T (S)')
ylabel('Y')
clc
output=[ArrayT',ArrayY',ArrayYTHEORY'];
save datfil output  -ascii
disp 'simulation finished'
