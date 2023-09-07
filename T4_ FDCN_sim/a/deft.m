function [f] = deft(q)
% Five-Point Endpoint Formula 4.7 burden & faires 9th ed
% a=[-25,48,-36,16,-3];q=length(x);
z=[3267,-29664,120008,-284256,435330,-448672,312984,-138528,29531]./5040;
a=wrev(z);
q1=length(a);d=q-q1+1;
c=[a,zeros(1,q-q1)];e=repmat(c,d,1);f1=0:d;g=e;
for i=1:d
f(i,:)=circshift(e(i,:),f1(i),2);
end
for i=d+1:q
  f(i,:)=[zeros(1,i-q1),z,zeros(1,q-i)];  
end
% f=f*x;
lk=1;
% g=flip(inv(f));
end
% b=[zeros(1,q-8),a;...
%     zeros(1,q-7),3,-16,36,-48,25,0,0;...
%     zeros(1,q-6),3,-16,36,-48,25,0;...
%     zeros(1,q-5),3,-16,36,-48,25];
