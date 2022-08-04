overlay = 1;
check = 'H8';
y = double(check(1)) - 64;
x = str2double(check(2:end));
Ind = name2(y,x);
autoROI1;
title(PlateCoordinate(Ind));