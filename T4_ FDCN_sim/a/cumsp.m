function CS=cumsp(x,y)
%CUMulative Simpson integration using SPLINE interpolation
%Call:
%      CS=cumspline(x,y)
%Input:
%		x = independent column
%		y= column-vect. or matrix of data columns
%Output:
%			
%	Vassili Pastushenko	 March 2006
%==============================
[ROWS,COLS]=size(y); 
N=2*ROWS-1;
X(1:2:N)=x;
X(2:2:N)=(x(1:end-1)+x(2:end))/2;
Y=interp1(x,y,X(:),'spline');
for i=1:3
 IND(:,i)=i:2:N+i-3;
end
FIL=[1 4 1]'/6;
  
CS=zeros(ROWS,COLS);
D=diff(x);

for i=1:COLS
    w=Y(:,i);
CS(2:ROWS,i)=(w(IND)*FIL).*D;
end

CS=cumsum(CS);


