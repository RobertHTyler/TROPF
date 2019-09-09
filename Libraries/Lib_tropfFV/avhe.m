function [dellav]=avhe(dell);
% Averages values with their eastern neighbor to give value at i+1/2,j,k
[Mm,Nn,Pp]=size(dell);
dellav=dell*0;
dellav(1:Mm-1,:,:)=(dell(1:Mm-1,:,:)+dell(2:Mm,:,:))/2;
dellav(Mm,:,:)=dellav(Mm-1,:,:);

