function [xtf,ytf]=predictpz(tf,xdum,ydum,x1dum,y1dum)
h=.01;
a=2.0926e7;
gm=1.4077e16;
t=0.;
x=xdum;
y=ydum;
x1=x1dum;
y1=y1dum;
while t<=(tf-.00001)
	xold=x;
	yold=y;
	x1old=x1;
	y1old=y1;
	step=1;
	flag=0;
	while step <=1
		if flag==1
			x=x+h*xd;
			y=y+h*yd;
			x1=x1+h*x1d;
			y1=y1+h*y1d;
			t=t+h;
			step=2;
		end
		tembot=(x^2+y^2)^1.5;
		x1d=-gm*x/tembot;
		y1d=-gm*y/tembot;
		xd=x1;
		yd=y1;
		flag=1;
	end;
	flag=0;
 	x=(xold+x)/2+.5*h*xd;
	y=(yold+y)/2+.5*h*yd;
	x1=(x1old+x1)/2+.5*h*x1d;
	y1=(y1old+y1)/2+.5*h*y1d;
end
xtf=x;
ytf=y;
