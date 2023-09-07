function [X,Y]=parfil(varargin)
%FILls midpoints into vectors/matrices by PARabolic interpolation
%as would do    Y=interp1(x,y,X,'parabolic')
%with X=sort([x;( x(1:end-1)+x(2:end) )/2];
%For matrices PARFIL works on columns. If only rows match length(x),
%then y will be transposed.
%For regular  data (equidistantly sampled, i.e. std(diff(x))=0), 
%Y is identical with INTERPOL(y,1,DEG)  in inner intervals,
%(completely identical for DEG=3), but INTERPOL is somewhat 
%faster than PARFIL
%Can have advantages if data are slightly noisy
%Call:
%       [X,Y]=parfil([x,]y[,DEG])
%       [X,Y]=parfil([DEG,][x,]y)
%       [X,Y]=parfil([x,][DEG,]y)
%Input:
%		x = strictly monotonic (i.e. distinct, all(diff(x)>0)=true)
%       vector of independent variable, length(x) = LENX >=3
%       default x=1:length(y);
%       y = vector of length LENX, or matrix, then any(size(y)==LENX)
%       If both x and y are vectors, x should precede y in input list
%
%       DEG = any scalar, but usually DEG=2 or DEG=3. Default DEG=2.
%       DEG is polynomial degree to fill the first and the last intervals.  
%
%       x,y,DEG may  appear in the list of input parameters in arbirary
%       sequence, but if x and y are vectors, then x should precede y.
%Output:
%       X,Y = interpolated data (i.e. containing midpoints)
%       X = column, Y = column if isvector(y).
%Subfunctions
%   [x,y,DEG]=getinput(varargin) GETINPUT checks input arguments
%   Y=parlagran(x3,y3,t) = vectorized parabolic and
%   Y=cublagran(x4,y4,t) = cubic elementary Lagrange interpolators
%   PARLAGRAN and CUBLAGRAN can be used for first and last interval(s)
%   PARFILDEMO: will called by PARFIL() (random x)
%Application area:
%  For very smooth data y, PARFIL(x,y,3) is close to SPLINE,
%  but it has great advantages over interp1(...'linear',or 'cubic',or
%  'spline') if the data are not densely sampled or have even a 
%   very small noise
% 
%   DEMO
%
%   x=(0:.1:20)';
%   x=2*(x+.09*(rand(size(x))-.5));
%   f=inline('cos(x).*exp(cos(x/sqrt(8)))');
%   y=[sin(x) f(x) exp(x/5)];       %Ideal data: SPLINE is the champ
%   y=y+.001*(rand(size(y))-.5); %with slightly noisy data: 
%                  %PARFIL can be the champ
%   [X,Y]=parfil(x,y);
%   YT=[sin(X) f(X) exp(X/5)]; %expected values
%   IND=3:199;
%   ERRPARF=std(Y(IND,:)-YT(IND,:));    %Error PARFIL
%   YPCHIP=interp1(x,y,X,'cubic'); %'cubic=pchip'
%   ERRPCHIP=std(YPCHIP(IND,:)-YT(IND,:))    %Error PCHIP
%   YINT1=interp1(x,y,X,'linear');  %Linear interpolation
%   ERRINT1=std(YINT1(IND,:)-YT(IND,:))  %Error linear   
%   YS=interp1(x,y,X,'spline'); 
%   ERRSPLINE=std(YS(IND,:)-YT(IND,:))  %Error spline
%   RES=[ERRPARF;ERRINT1; ERRPCHIP; ERRSPLINE];
%    
%   METHODS={'PARFIL    ';'LINEAR         ';'PCHIP_CUB  ';'SPLINE         '};
%   disp('                  1e6* ERRORS')
%   for i=1:4
%   disp([METHODS{i},'  ',sprintf('%15.2f',1e6*RES(i,:))])
%   end

%	Vassili Pastushenko	 5-th March 2006
%==============================
ARGS=varargin;
if isempty(ARGS)
    parfildemo
    return
end
[x,y,DEG]=getinput(ARGS);

 [N,C]=size(y);
 Nm1=N-1;
 Nm2=N-2;
 Nm3=N-3;
 N2m1=N*2-1; 
 
 X(1:2:N2m1,1)=x;
 X(2:2:N2m1,1)=(x(1:Nm1)+x(2:N))/2;

 %Indices
 indNm2=1:Nm2;
 indNm3=1:Nm3;
 ods=1:2:N2m1;
 
 ind2Nm1=2:Nm1;
 ind2Nm2=2:Nm2;
 evs=2:2:N2m1;
 
 dx=diff(x);
 DEL=dx(1:Nm2)+dx(2:Nm1);
 T(:,1)= dx(indNm2)./DEL; %interp [0 dy1 dy2] in [0 T 1]
 MUL=.25./(T-1);
 FILL=[MUL.*(T-2),MUL.*T.^2];
 MUR=.25*(T+1);
 FILR=[MUR./T, MUR];
 Y=zeros(N2m1,C);
 
 %Processing y
 
 for col=1:C
 yw=y(:,col);
 dy=diff(yw);
 AD=yw(indNm2);
 % Insertions
 DY=[dy(indNm2),dy(indNm2)+dy(ind2Nm1)];
 Lefty=sum(FILL.*DY,2)+AD;
 Righty=sum(FILR.*DY,2)+AD;
 
 WLEF=X(4:2:N2m1-3)-X(1:2:N2m1-6);
 WRIG= X(7:2:N2m1)-X(4:2:N2m1-3);
 
 WLEF=WLEF./(WLEF+WRIG);
 WRIG=1-WLEF;
 if N2m1>5
 InsY=[Lefty(1);  (Lefty(ind2Nm2).*WLEF+WRIG.*Righty(indNm3)); Righty(Nm2)];
 else
 InsY=[Lefty(1); Righty(Nm2)];
 end
 
 %Combine
 Y(ods,col)=yw;
 Y(evs,col)=InsY; 
 end
 if N2m1==5||all(DEG~=(2:3))
     return
 end
 
 switch DEG
     case 2
 inn=[1,3:4];
 Y(2,:)=parlagran(X(inn),Y(inn,:),X(2));
 inn=N2m1-[3 2 0];
 rep=N2m1-1;
 Y(rep,:)=parlagran(X(inn),Y(inn,:),X(rep));
     case 3
 inn=[1,3:5];
 Y(2,:)=cublagran(X(inn),Y(inn,:),X(2));
 inn=N2m1-[4 3 2 0];
 rep=N2m1-1;
 Y(rep,:)=cublagran(X(inn),Y(inn,:),X(rep));
 end 
 
 function Y=parlagran(x,y,t)
%Y=interp1(x,y,t,'parabolic'), x(1)<t<x(3) scalar, size(x)=[3 1],
%size(y)= [3 C], C>=1
%  size(Y)=[1 C]
      yy=y-repmat(y(1,:),3,1);
      t=t-x(1);
      x=x-x(1);
      x2=x(2);x3=x(3);
     %Lagrange parabolic origin   
 Y=t*(t-x3)/(x2*(x2-x3))*yy(2,:) + t*(t-x2)/(x3*(x3-x2)) * yy(3,:) +y(1,:);
  %Y=(yy(2,:).*(x(3)-t)*(t/x(2))+yy(3,:).*(t-x(2))*(t/x(3)))/(x(3)-x(2) )+y(1,:); %parlagran     
 
 
 function Y=cublagran(x,y,t)
%Y=interp1(x,y,t,'parabolic'), x(1)<t<x(4) scalar, size(x)=[4 1],
%size(y)= [4 C], C>=1
%  size(Y)=[1 C]
      yy=y-repmat(y(1,:),4,1);
      t=t-x(1);
      x=x-x(1);
      x2=x(2);x3=x(3);x4=x(4);
   %Lagrange cubic origin   
   Y=t*(t-x3)*(t-x4)/(x2*(x2-x3)*(x2-x4)) * yy(2,:)+...
     t*(t-x2)*(t-x4)/(x3*(x3-x2)*(x3-x4)) * yy(3,:)+...
     t*(t-x2)*(t-x3)/(x4*(x4-x2)*(x4-x3)) * yy(4,:) +y(1,:);
  
 function [x,y,DEG]=getinput(ARGS)
 LARG=length(ARGS);
 switch LARG
     
     case 1   %Only y present
         y=ARGS{1};
     case 2
         for i=1:2
            LENS(i)=length(ARGS{i});
         end
         NU=find(LENS==1);
         
         if ~isempty(NU)
             NU=NU(1);
             DEG=ARGS{NU};
             ARGS(NU)=[];
             y=ARGS{1};
         else
             x=ARGS{1};
             y=ARGS{2};
             if ~isvector(x)
                 t=y;
                 y=x;
                 x=t;
             end
         end
     case 3
         for i=1:3
            LENS(i)=length(ARGS{i});
         end
         NU=find(LENS==1);
         NU=NU(1);
         if isempty(NU)
             o_o;
             pause(.2)
             error('DEG should be scalar')
         end
             DEG=ARGS{NU};
             ARGS(NU)=[];
             x=ARGS{1};
             y=ARGS{2};
             if ~isvector(x)
                 t=y;
                 y=x;
                 x=t;
             end
     otherwise
         o_o,pause(.2)
         error('Too many input parameters')
 end
     
 if isvector(y)
     y=y(:);
 end
 
 if size(y,1)<3
     y=y.';
     if size(y,1)<3
         o_o
         error('Too short data')
     end
 end
 
 
 if ~exist('x')
     x=(1:size(y,1));
 end
 x=x(:);
 LENX=length(x);
 
 if LENX<3
     checkpoint
     error('Too short x')
 end
 [R,C]=size(y);
 
 if LENX~=R
     if LENX==C
         y=y.';
     else
         o_o;
         pause(.2)
         error('size(x) mismatches size(y)')
     end
 end     
 
 if ~exist('DEG')
     DEG=2;
 end

 
 function checkpoint
    o_o;
    display('Check data');
    pause(.2)
function parfildemo
%
%Call:
%     
%Input:
%		
%		
%Output:
%			
%	Vassili Pastushenko	 March 2006
%==============================
  
    x=(0:.1:10)';
    x=2*(x+.08*(rand(size(x))-.5));
    f=inline('cos(x).*exp(cos(x/sqrt(8)))');
    y=[sin(x) f(x) exp(x/10)];       %Ideal data: SPLINE is the champ
    y=y+.005*(rand(size(y))-.5); %with slightly noisy data: 
                   %PARFIL can be the champ
    [X,Y]=parfil(x,y);
    YT=[sin(X) f(X) exp(X/10)]; %expected values
    IND=4:198;
    ERRPARF=std(Y(IND,:)-YT(IND,:));    %Error PARFIL
    YPCHIP=interp1(x,y,X,'cubic'); %'cubic=pchip'
    ERRPCHIP=std(YPCHIP(IND,:)-YT(IND,:))    %Error PCHIP
    YINT1=interp1(x,y,X,'linear');  %Linear interpolation
    ERRINT1=std(YINT1(IND,:)-YT(IND,:))  %Error linear   
    YS=interp1(x,y,X,'spline'); 
    ERRSPLINE=std(YS(IND,:)-YT(IND,:))  %Error spline
    RES=[ERRPARF;ERRINT1; ERRPCHIP; ERRSPLINE];
     
    METHODS={'PARFIL    ';'LINEAR         ';'PCHIP_CUB  ';'SPLINE         '};
    disp('                  1e6* ERRORS')
    for i=1:4
    disp([METHODS{i},'  ',sprintf('%15.2f',1e6*RES(i,:))])
    end
    subplot(2,2,1)
    cla
    setcol
     plot(X,YT(:,1:2),'.-');
     hold on
    plot(x,y(:,1:2),'k.');
    
    set(gca,'fontsize',15)
    title('PARFIL data: black points')
    xlabel('irregular x')
    ylabel('sin(x),   cos(x).*exp(cos(x/sqrt(8))')
    hold off
    axis tight
    
    subplot(2,2,4)
    cla
    setcol
    hold on
    plot(X(IND),1e3*(Y(IND,1:2)-YT(IND,1:2)),'.',X(IND),1e3*(YPCHIP(IND,1:2)-YT(IND,1:2)),'.')
    set(gca,'fontsize',15)
    plot(X(IND),1e3*(Y(IND,1:2)-YT(IND,1:2)),'.')
    legend('PARFIL','PARFIL','PCHIP','PCHIP')
    plot(X(IND),1e3*(Y(IND,1:2)-YT(IND,1:2)),':',X(IND),1e3*(YPCHIP(IND,1:2)-YT(IND,1:2)),':')
    title('1000 Errors')
    hold off
    set(gca,'fontsize',15)
    axis tight
    shg
    function setcol
    M=get(gca,'colororder');
    M(1:4,:)=[1 0 0
        0 .8 0
        0 0 1
        .7 0 .7];
   
    set(gca,'colororder',M);
    hold on
