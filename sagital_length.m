function [ length ] = sagital_length()
%-------------------------------------------------------------------------
% This function estimates the aortic arch length from sagittal plane image
% Output
% length: estimated aortic arch length
%-------------------------------------------------------------------------

[filename, pathname] = uigetfile('*.*');

filename = fullfile (pathname, filename);
info = dicominfo(filename);
Y = dicomread(info);

coordx = zeros(10,1);
coordy = zeros(10,1);
xs = [];
ys = [];
xold=0;
yold=0;
imshow(Y, 'DisplayRange', []);hold on 

for k = 1:10
    [coordx(k,1), coordy(k,1)] = ginput(1);  
    xs = [xs;coordx(k,1)];
    ys = [ys;coordy(k,1)];
    if xold;
        plot([xold coordx(k,1)],[yold coordy(k,1)],'r.-');
    else
        plot(coordx, coordy,'go');
    end 
    xold=coordx(k,1);
    yold=coordy(k,1);
    %coordinates(k,:) = [coordx(k,1), coordy(k,1)];
    %plot(coordinates(:,1), coordinates(:,2), '.-');
end

ps = info.PixelSpacing(1); % pixel spacing (from DICOM)
length = 0; % initializes the length of the aorta

for k = 1 : 9
    segment = sqrt(((coordx(k,1)-coordx(k+1,1))*ps)^2 + ((coordy(k,1)-coordy(k+1,1))*ps)^2);
    length = length + segment;
end

end

