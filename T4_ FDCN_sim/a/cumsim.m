function Y=cumsim(varargin)
% CUMulative SIMpson integration 
% of vectors or matrices (usually columns)
%
%Call:
%       Y=cumsim(x,y,DEG)
%Input:
%		y=data(vector or matrix)
%		x = vector of independent variable (all(diff(x)>0)=true, length(x)>=3) 
%       by default x=1:eval('S=size(y), S=S(S>2)')
%       if set, any(length(x)==size(y))=true
%
%       DEG = (def. 3) argument for PARFIL(x,y,DEG) 
%       used for interpolation 
%Output:
%		Approximation to cumulative Simpson integral, ordered in columns
%       Y = integral(y(x),x(1),x)

%	Vassili Pastushenko	 23.02.2006
%==============================
%Set default DEG = 3 
ARGS=varargin;
 DEG=[];
 for i=1:length(ARGS)
     if length(ARGS{i})==1
         DEG=ARGS{i};
     end
 end
 if isempty(DEG)
     ARGS{end+1}=3;
 end
     [x,y,DEG]=getinput(ARGS);
     LENDAT=size(y,1);
     dx=diff(x);
 %Integration
 
 [XX,YY]=parfil(x,y,DEG); %Parabolic insertion midpoints
 [R,C]=size(YY);
 for i=1:3
 IND(:,i)=i:2:R+i-3;
 end
 
 
 FIL=[1 4 1]'/6;
 Y=zeros(LENDAT,C);
 IN=2:LENDAT;
 
 for col=1:C 
 yw=YY(:,col);
 Y(IN,col)=(yw(IND)*FIL).*dx;
 end
 Y=cumsum(Y);
 
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
