function ER=compcum
%comparison of different cumulative Simpson sums
%for regular or irregular exponential data
%Call:
%        ER=compcum(reg)
%Input:
%		if set, regular  abscissa (equally spaced)
%		otherwise irregular spacing
%Output:
%		ER= [absolute errors; relative errors]
%       methods: [CUMSIM, CUMSIMPSON, CUMSIMPSUM]
%
%	Vassili Pastushenko	 March 2006
%====================================================
 
NP=100; %number of points
x=4*pi*sort(rand(NP,1)); %irregular


com(x,1); %DEMO error increments
com(x,2); %DEMO error cumsum
figure(gcf)

function com(x,sif)
MUL=10^(6-(sif>1));
subplot(1,2,sif);

y=sin(x);  %simulated data
yt=[cos(x(1))-cos(x)];  %expected answer: yt=integral(y,0,x)
c=cumsim(x,y);              %CUMSIM
z=cumsimpson(x,y);         %CUMSIMPSON, cf. File Exchange
zz=cumspline(x,y);          %CUMSPLINE
Z=[c z zz];
ZT=repmat(yt,1,3);
ERR=MUL*(Z-ZT);
if sif<2
ERR=diff(ERR);
x=parfil(x,x);
x=x(2:2:end);
end

err=sqrt(mean(ERR.^2));
cla
setcol
plot(x,ERR,'.','markersize',10);
axis tight
legend('CUMSIM','CUMSIMPSON','CUMSPLINE','location','best');
plot(x,ERR,':');
hold off
set(gca,'fontsize',15)
if sif<2
ylabel('diff(CS(x,sin(x)))-diff(cos(x))')
else
    ylabel('CS(x,sin(x))-cos(x)+1')
end

xlabel('x')

switch sif
    case 1
TEX=['std(increment errors) '];
    case 2
        TEX=['std(errors) '];
end
    
TEX={TEX;['[',sprintf('%8.1f',err),']/',num2str(MUL)]};
title(TEX);

function setcol
M=get(gca,'colororder');
M(1:3,:)=[1 0 0;0 .9 0;0 0 1];
set(gca,'colororder',M);
hold on