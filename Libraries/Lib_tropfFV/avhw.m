function [dellav]=avhw(dell);
% Averages values with their western neighbor to give value at i-1/2,j,k
[Mm,Nn,Pp]=size(dell);
dellav=dell*0;
dellav(2:Mm,:,:)=(dell(1:Mm-1,:,:)+dell(2:Mm,:,:))/2;
dellav(1,:,:)=dellav(2,:,:);

