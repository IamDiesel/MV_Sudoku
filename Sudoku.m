%% clear all
imtool close all;
close all;
clear;

%% read image
%I = imread('../sudoku3.jpg'); %read rgb image (photo)
%I = imread('../sudoku1.tif'); %read rgb image (cut photo)
I = imread('../sudoku1.png'); %read rgb image (wikipedia)
imtool(I);

%% convert to gray image
IG = rgb2gray(I);
imtool(IG);

%% calculate gray threshold
figure;
imhist(IG);
[counts,binLocations]=imhist(IG);
countsFilt = medfilt2(counts,[50 1]);%smooth histogramm
LMax = imregionalmax(countsFilt);%find maxima
LMin = imregionalmin(countsFilt);%find minima
LMax = LMax(2:length(LMax)) - LMax(1:length(LMax)-1);%derivate of maxiums
LMin = LMin(2:length(LMin)) - LMin(1:length(LMin)-1);%derivate of minimums
LMax1 = find(LMax,2);%first maximum
LMin1 = find(LMin(LMax1(2):length(LMin)),2);%first minimum after first maximum
if length(LMin1) > 1
    [countMin, binLocationMin] = min(counts(LMax1(2)+LMin1(1):LMax1(2)+LMin1(2)));
else
    [countMin, binLocationMin] = min(counts(LMax1(2)+LMin1(1): length(counts)));
end
binLocationMin = binLocationMin + LMax1(2) + LMin1(1) - 2;%-2 because of index displacement
grayThreshold = double(binLocationMin/255);

%% create BW image
BW = imcomplement(im2bw(IG,grayThreshold)); %gray -> bw + invert
imtool(BW);
imageArea = numel(BW);
imageSize = size(BW);

%% remove noise

BW = bwareaopen(BW,ceil(imageArea/5000));
imtool(BW);

%% clear border
%BW = imclearborder(BW);
%imtool(BW);

%% %% find the largest object (box)
hold on;
R = regionprops(BW,'Area','BoundingBox','PixelList');
NR = numel(R);

maxArea = 0;
for k = 1:NR
    A(k) = prod(R(k).BoundingBox(3:4)); %area of BoundingBox
    if R(k).Area > maxArea
        maxArea = R(k).Area;
        kmax = k;
    end
end


BBmax = R(kmax).BoundingBox; % biggest bounding box
sumXY = sum(R(kmax).PixelList,2);% sum of x and y pixels
diffXY = diff(R(kmax).PixelList,[],2);% diff of x and y pixels

[m,dUL] = min(sumXY); % upper left corner: sum of x and y min
[m,dLR] = max(sumXY); % lower right corner: sum of x and y max
[m,dLL] = min(diffXY); % lower left corner: diff of x and y min
[m,dUR] = max(diffXY); % upper right corner: diff of x and y max

pts = R(kmax).PixelList([dUL dLL dLR dUR dUL],:);
figure;
imshow(I);
hold on;
h_pts = plot(pts(:,1),pts(:,2),'m','linewidth',3);

XYLIMS = [BBmax(1) + [0 BBmax(3)] BBmax(2) + [0 BBmax(4)]];

%% identify objects inside the box
A_tmin = imageArea/5000; % Bounds for the digit pixel area
A_tmax = imageArea/150;
digitbox_minarea = imageArea/2500; % Bounds for the digit bounding box area
digitbox_maxarea = imageArea/75;%9x9 fields: imageSize/91

kgood = zeros(1,NR);
Pnew = zeros(NR,2);
for k = 1:NR
	if R(k).Area < A_tmax && A(k) > digitbox_minarea && A(k) < digitbox_maxarea ...
        && R(k).BoundingBox(3) < imageSize(1)/8 && R(k).BoundingBox(4) < imageSize(2)/8 ...
        && R(k).BoundingBox(3) > 1 && R(k).BoundingBox(4) > 1
                
            Pnew(k,:) = [R(k).BoundingBox(1)+R(k).BoundingBox(3)/2 R(k).BoundingBox(2)+R(k).BoundingBox(4)/2];       
            if inpolygon(Pnew(k,1),Pnew(k,2),pts(:,1),pts(:,2))
                h_digitcircles(k) = plot(Pnew(k,1),Pnew(k,2),'ro','markersize',24);
            end      
    end
end

%% Draw the grid based on the corners
T = cp2tform(pts(1:4,:),0.5 + [0 0; 9 0; 9 9; 0 9],'projective');
for n = 0.5 + 0:9, [x,y] = tforminv(T,[n n],[0.5 9.5]); plot(x,y,'g'); end
for n = 0.5 + 0:9, [x,y] = tforminv(T,[0.5 9.5],[n n]); plot(x,y,'g'); end

%% Only keep elements in the boxes
T = cp2tform(pts(1:4,:),[0.5 0.5; 9.5 0.5; 9.5 9.5; 0.5 9.5],'projective');
Plocal = (tformfwd(T,Pnew));
Plocal = round(2*Plocal)/2;

del = find(sum(Plocal - floor(Plocal) > 0 |  Plocal < 1 | Plocal > 9,2)) ;
Pnew(del,:) = [];

%delete(nonzeros(h_digitcircles(del)));
delete(h_digitcircles(del));

%% Show the coordinate transforms
figure;
T = cp2tform(pts(1:4,:),500*[0 0; 1 0; 1 1; 0 1],'projective');
IT = imtransform(double(BW),T);
imshow(IT);

%% Show the template data
figure;
load TEMPLATEDATA
for n = 1:9
    subplot(3,3,n),imagesc(NT{n});
end
colormap gray;

%% Calculate the Solution


Plocal = identifynumbers_fun(pts,Pnew,NT,BW);
M = zeros(9);
for k = 1:size(Plocal,1)
    M(Plocal(k,2),Plocal(k,1)) = Plocal(k,3);
end
M_sol = drawgraph(M);



% remove big white areas (percentage of white pixels in rectangular space of object high)
% CC = bwconncomp(BW);
% numPixels = cellfun(@numel,CC.PixelIdxList);
% imagesize = CC.ImageSize(1) * CC.ImageSize(2);
% for i = 1:length(numPixels)
%     x = floor(CC.PixelIdxList{i}/CC.ImageSize(1));
%     y = mod(CC.PixelIdxList{i},CC.ImageSize(1));
%     if numPixels(i)/((max(x)-min(x))*(max(y)-min(y))) > 0.167 && ((max(x)-min(x))*(max(y)-min(y)))/imagesize > 0.1
%         BW(CC.PixelIdxList{i}) = 0;
%     end
% end
%  imtool(BW);
% 
% BWT=bwmorph(BW, 'thin', inf);
% imtool(BWT);
% 
% [H, T, R] = hough(BWT);
% figure;
% imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,...
%       'InitialMagnification','fit');
% title('Hough Transform');
% xlabel('\theta'), ylabel('\rho');
% axis on, axis normal, hold on;
% colormap(hot);
% NHS=floor(((size(H)/20)-1)/2)*2+1; %NHoodSize Parameter
% P  = houghpeaks(H,99,'threshold',ceil(0.3*max(H(:))),'NHoodSize',NHS);
% x = T(P(:,2)); y = R(P(:,1));
% plot(x,y,'s','color','white');
% % Find lines and plot them
% FG = length(BW); %FillGap Parameter
% lines = houghlines(BWT,T,R,P,'FillGap',FG,'MinLength',7);
% % figure, imshow(RGB), hold on
% for k = 1:length(lines)
%    xy = [lines(k).point1; lines(k).point2];
% %    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% 
%    % Plot beginnings and ends of lines
% %    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
% %    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% 
% end
% 
% figure, imshow(I), hold on
% for i = 1:length(lines)
%    line(i,1) = lines(i).point1(1);
%    line(i,2) = lines(i).point1(2);
%    line(i,3) = lines(i).point2(1) - lines(i).point1(1);
%    line(i,4) = lines(i).point2(2) - lines(i).point1(2);
%    plot([line(i,1) line(i,1)+line(i,3)],[line(i,2) line(i,2)+line(i,4)],'LineWidth',2,'Color','green')
% end
% n=1;
% for i = 1:length(lines)
%     for j = i+1:length(lines)
%        r1  = line(i,3)/line(j,3);
%        r2 = line(i,4)/line(j,4);
%        if (r1 ~= r2 && i ~=j)
%            A = [line(i,3) -line(j,3);
%                 line(i,4) -line(j,4)]; 
%            b = [line(j,1)-line(i,1);
%                line(j,2)-line(i,2)];
%            temp = A\b;
%            if ((abs(temp(1)) ~= Inf) || (abs(temp(1)) ~= Inf))
%                 L(:,n) = [i;
%                           j;
%                           line(i,1)+temp(1)*line(i,3);
%                           line(i,2)+temp(1)*line(i,4)];
%                 plot(L(3,n),L(4,n),'x','LineWidth',2,'Color','red');
%                 n = n+1;
%            end
%        end
%     end
% end