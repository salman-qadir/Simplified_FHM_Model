 function sz=sq(mn,initial,cr)
 
  w=par;   sz=initial(w,cr); 
  
for ij=1:(length(sz.tsp)-1) 
    
     sz=mn(ij,w,sz);    
     if rem(ij,1000/w.t)==0
        fprintf('\rij\t=\t%5.1f',ij); 
        fprintf('\rVolt\t=\t%5.2f',sz.v(ij,3));
     end     
     if (sz.v(ij,3)<2.5)||(min(sz.ce(ij,:))<10)||sz.SOC(ij+1,2)<.1||sz.SOC(ij+1,1)<.1          
        fprintf('\rij\t=\t%1f',ij); 
        fprintf('\rVolt\t=\t%5.2f',sz.v(ij,3));
        fprintf('\rSOCN\t=\t%5.2f',sz.SOC(ij+1,1));
        fprintf('\rSOCP\t=\t%5.2f',sz.SOC(ij+1,2));
        fprintf('\rCE\t=\t%5.2f\r',min(sz.ce(ij,:)));
        break;
     end    
end
sz.t=sz.tsp(1:ij);
end
 
 
%  matname=fullfile('C:\Users\exx241\Dropbox\Tasz\T4\shm','sz.mat');
%  delete sz.mat
%  save(matname,'sz');
%  clearvars;
%  load sz;
% m = readmatrix('set.csv');v1=m(:,1);v2=m(:,2);
% if length(sz.v(:,3))>length(v1)
% vb=interp(v1,round(length(sz.v(:,3))/length(v1)));
% va=interp(v2,round(length(sz.v(:,3))/length(v1)));
% v3=sz.v(:,3);
% else
% v3=interp(sz.v(:,3),ceil(length(sz.v(:,3))\length(v1)));
% vb=v1;va=v2;
% end
% plot(vb,va);hold on;plot(vb,v3(1:length(vb)));grid on;