function [dellav]=avhs(dell);
% Averages values with their southern neighbor to give value at i,j-1/2,k
[Mm,Nn,Pp]=size(dell);
dellav=dell*0;
dellav(:,2:Nn,:)=(dell(:,1:Nn-1,:)+dell(:,2:Nn,:))/2;
dellav(:,1,:)=dellav(:,2,:);

