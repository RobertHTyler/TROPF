function [dellav]=avhn(dell);
% Averages values with their northern neighbor to give value at i,j+1/2,k
[Mm,Nn,Pp]=size(dell);
dellav=dell*0;
dellav(:,1:Nn-1,:)=(dell(:,1:Nn-1,:)+dell(:,2:Nn,:))/2;
dellav(:,Nn,:)=dellav(:,Nn-1,:);

