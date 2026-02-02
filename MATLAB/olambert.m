% olambert.m

function [vrx,vry]=olambert(xic,yic,tfdes,xf,yf,xlongm,xlongt)

global a gm % In version 3 comment this out!
global count ArrayICOUNT ArrayBGAM ArrayVRX ArrayVRY ArrayTF % for output array (if reqd)

% Initialise outputs (if reqd)
count=0;
vrx=0;
vry=0;
tf=0;

%A
ric=sqrt(xic^2+yic^2);
rf=sqrt(xf^2+yf^2);
cphi=(xic*xf+yic*yf)/(ric*rf);
phi=acos(cphi);
r0=ric;
degrad=360./(2.*pi);


% Initialise while loop
SecondTimeThrough=0;

% Program executes this loop twice
while SecondTimeThrough <= 1
   	% Initialise for loop
	if SecondTimeThrough == 0
		start=-90;
		step=.1;
		stop=+90;
	else
		start=gamdegnew;
		step=.0001;
		stop=gamdegfin;
	end;

	% Main body of program
	for gamdeg=start:step:stop

		%B

		gam=gamdeg/degrad;
		top=gm*(1-cos(phi));
		temp=r0*cos(gam)/rf-cos(phi+gam);
		bot=r0*cos(gam)*temp;


		if ~(top<0. | bot<0.)

			%C

			v=sqrt(top/bot);
			if xlongt>xlongm
		  		vrx=v*cos(pi/2 -gam+xlongm);
		  		vry=v*sin(pi/2 -gam+xlongm);
			else
		  		vrx=v*cos(-pi/2 +gam+xlongm);
		  		vry=v*sin(-pi/2 +gam+xlongm);
			end
			xlam=r0*v*v/gm;
			top1=tan(gam)*(1-cos(phi))+(1-xlam)*sin(phi);
			bot1p=(1-cos(phi))/(xlam*cos(gam)*cos(gam));
			bot1=(2-xlam)*(bot1p+cos(gam+phi)/cos(gam));
			top2=2*cos(gam);

			if ~((2/xlam-1) < 0.)

				%D

				bot2=xlam*((2/xlam-1)^1.5);
				top3=sqrt(2/xlam-1);
				bot3=cos(gam)/tan(phi/2)-sin(gam);
				temp=(top2/bot2)*atan2(top3,bot3);
				tf=r0*(top1/bot1+temp)/(v*cos(gam));

				if (tf > tfdes)

					break % out of the for loop
				end; % condition #3

				% output arrays (if reqd)

   				count=count+1;
				ArrayICOUNT(count)=count;
  				ArrayBGAM(count)=57.3*gam;
  				ArrayVRX(count)=vrx;
  				ArrayVRY(count)=vry;
  				ArrayTF(count)=tf;


			end; % condition #2
		end; % condition #1

	end; % for loop


%E

gamdegnew=gamdeg-.15;
gamdegfin=gamdeg+1.;

SecondTimeThrough=SecondTimeThrough+1;
end % while loop
plot(count,ArrayTF)
