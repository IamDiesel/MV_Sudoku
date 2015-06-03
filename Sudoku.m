%% clear all
imtool close all;
close all;
clear;

%% read image
I = imread('../sudoku6.jpg'); %read rgb image (photo)
% I = imread('../sudoku1.tif'); %read rgb image (cut photo)
% I = imread('../sudoku1.png'); %read rgb image (wikipedia)
% I = imread('../sample.bmp');
% imtool(I); % original image

%% convert to gray image

if numel(size(I)) == 3
    IG = rgb2gray(I);

else
    IG = I;
end

% imtool(IG); % grayscale image
imageArea = numel(IG); % number of pixels
imageSize = size(IG); % width and height

%% calculate gray threshold
% figure;
% imhist(IG);
% [counts,binLocations]=imhist(IG);
% countsFilt = medfilt2(counts,[50 1]);%smooth histogramm
% LMax = imregionalmax(countsFilt);%find maxima
% LMin = imregionalmin(countsFilt);%find minima
% LMax = LMax(2:length(LMax)) - LMax(1:length(LMax)-1);%derivate of maxiums
% LMin = LMin(2:length(LMin)) - LMin(1:length(LMin)-1);%derivate of minimums
% LMax1 = find(LMax,2);%first maximum
% LMin1 = find(LMin(LMax1(2):length(LMin)),2);%first minimum after first maximum
% if length(LMin1) > 1
%     [countMin, binLocationMin] = min(counts(LMax1(2)+LMin1(1):LMax1(2)+LMin1(2)));
% else
%     [countMin, binLocationMin] = min(counts(LMax1(2)+LMin1(1): length(counts)));
% end
% binLocationMin = binLocationMin + LMax1(2) + LMin1(1) - 2;%-2 because of index displacement
% grayThreshold = double(binLocationMin/255);

grayThreshold = median(double(IG(:)))/1.3/255; % grayThreshold

%% create BW image

BW = imcomplement(im2bw(IG,grayThreshold)); %gray -> bw + invert
% imtool(BW); % black and white image

%% remove noise

BW = bwareaopen(BW,ceil(imageArea/5000)); % opening
% imtool(BW); % black and white image + opening

%% clear border
BW = imclearborder(BW);
% imtool(BW); % black and white image + opening + removing white objects
% frpm border (e.g. shadows in the image)

%% %% find the largest object (box)

R = regionprops(BW,'Area','BoundingBox','PixelList');
NR = numel(R);
kamx = 0;
maxArea = 0;
for k = 1:NR
    if R(k).Area > maxArea
        maxArea = R(k).Area;
        kmax = k;
    end
end

sumXY = sum(R(kmax).PixelList,2);% sum of x and y pixels
diffXY = diff(R(kmax).PixelList,[],2);% diff of x and y pixels

[m,dUL] = min(sumXY); % upper left corner: sum of x and y min
[m,dLR] = max(sumXY); % lower right corner: sum of x and y max
[m,dLL] = min(diffXY); % lower left corner: diff of x and y min
[m,dUR] = max(diffXY); % upper right corner: diff of x and y max

pts = R(kmax).PixelList([dUL dLL dLR dUR dUL],:);

%% transform image and cut sudoku

T = fitgeotrans(pts(1:4,:),900*[0 0; 1 0; 1 1; 0 1],'projective');
IT = logical(imwarp(double(BW),T));
R = regionprops(IT,'Area','BoundingBox','PixelList');
NR = numel(R);
kmax = 0;
maxArea = 0;
for k = 1:NR
    A(k) = prod(R(k).BoundingBox(3:4)); %area of BoundingBox
    if R(k).Area > maxArea
        maxArea = R(k).Area;
        kmax = k;
    end
end
BBmax = R(kmax).BoundingBox;
IT = IT(ceil(BBmax(2)):ceil(BBmax(2)+BBmax(4)),ceil(BBmax(1)):ceil(BBmax(1)+BBmax(3)));
imageArea = numel(IT); % number of pixels
imageSize = size(IT); % width and height


R = regionprops(IT,'Area','BoundingBox','PixelList');
NR = numel(R);
kmax = 0;
maxArea = 0;
for k = 1:NR
    A(k) = prod(R(k).BoundingBox(3:4)); %area of BoundingBox
    if R(k).Area > maxArea
        maxArea = R(k).Area;
        kmax = k;
    end
end
BBmax = R(kmax).BoundingBox;

XYlimits = [floor(BBmax(1)) + [1 BBmax(3)] floor(BBmax(2)) + [1 BBmax(4)]];
pts = [ XYlimits(3) XYlimits(1); XYlimits(3) XYlimits(2); XYlimits(4) XYlimits(2); XYlimits(4) XYlimits(1); XYlimits(3) XYlimits(1) ];
figure;
imshow(IT);
hold on;
h_pts = plot(pts(:,1),pts(:,2),'m','linewidth',3);

%% Draw the grid based on the corners

for n = linspace(XYlimits(1), XYlimits(2), 10)
    if n ~= XYlimits(1) && n ~= XYlimits(2)
        x = [ n n ];
        y = [ XYlimits(3) XYlimits(4) ];
        plot(y, x, 'g', 'linewidth', 2); 
    end
end
for n = linspace(XYlimits(3), XYlimits(4), 10)
    if n ~= XYlimits(1) && n ~= XYlimits(2)
        x= [ XYlimits(1) XYlimits(2) ]; 
        y = [ n n ]; 
        plot(y, x, 'g', 'linewidth', 2); 
    end
end

%% identify objects inside the box

A_tmin = imageArea/3000; % Bounds for the digit pixel area
A_tmax = imageArea/90;
digitbox_minarea = imageArea/2500; % Bounds for the digit bounding box area
digitbox_maxarea = imageArea/75;%9x9 fields: imageSize/81

kgood = zeros(1,NR);
Pnew = zeros(NR,2);
for k = 1:NR
	if R(k).Area < A_tmax && A(k) > digitbox_minarea && A(k) < digitbox_maxarea ...
        && R(k).BoundingBox(3) < imageSize(1)/8 && R(k).BoundingBox(4) < imageSize(2)/8 ...
        && R(k).BoundingBox(3) > 1 && R(k).BoundingBox(4) > 1
            
            kgood(k) = 1;
            Pnew(k,:) = [R(k).BoundingBox(1)+R(k).BoundingBox(3)/2 R(k).BoundingBox(2)+R(k).BoundingBox(4)/2];       
            if inpolygon(Pnew(k,1),Pnew(k,2),pts(:,1),pts(:,2))
                h_digitcircles(k) = plot(Pnew(k,1),Pnew(k,2),'ro','markersize',max(R(k).BoundingBox(3:4)));
            end      
    end
end

%% identify numbers






