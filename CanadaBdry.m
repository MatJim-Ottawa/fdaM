
%  load the lat (deg,min) and lon (deg,min) from csv file

load CanadaLatLon.csv
size(CanadaLatLon)

%  convert minutes to fractions

CanadaLatLon(:,2) = CanadaLatLon(:,2)./60;
CanadaLatLon(:,4) = CanadaLatLon(:,4)./60;

%  set up CanadaBdry, first  column is minus fractionated longitude,
%                     second column is fractionated latitude

CanadaBdry = zeros(62,2);
CanadaBdry(:,1) = -(CanadaLatLon(:,3)+CanadaLatLon(:,4));
CanadaBdry(:,2) =   CanadaLatLon(:,1)+CanadaLatLon(:,2);

%  save boundary points as a .txt file

save CanadaBdry.txt CanadaBdry -ascii

%  set up edge points

load CanadaBdry.txt
tmp = zeros(62,2);
for i=1:62
    tmp(i,1) = i;
    tmp(i,2) = i+1;
end
tmp(62,2) = 1;
save CanadaEdge.txt tmp -ascii
