% Preallocation
clear
Z=zeros(1,100);
count=0;
Z1=0;
for I=1:100
    SUM=0;
    for J=1:12
        RAND=rand(1);
        SUM=SUM+RAND;
    end;
    X=SUM-6;
    Z(I)=X;
    Z1=Z(I)+Z1;
    XMEAN=Z1/I;
end;
SIGMA=0;
Z1=0;
for I=1:100
    Z1=(Z(I)-XMEAN)^2+Z1;
    if I == 1
        SIGMA=0;
    else
        SIGMA=sqrt(Z1/(I-1));
    end;
    count=count+1;
    ArrayI(count)=I;
    ArraySIGMA(count)=SIGMA;    
end;
figure
plot(ArrayI,ArraySIGMA),grid
title('Sampled Standard Deviation')
xlabel('Number of Samples')
ylabel('Calculated Standard Deviation')
clc
disp '*** Simulation Complete'
output=[ArrayI',ArraySIGMA'];
save datfil.txt output  -ascii
disp 'simulation finished'
