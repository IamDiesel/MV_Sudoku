close all;
RGB = imread('../sudoku3.jpg'); %read rgb image (photo)
RGB = imread('../sudoku1.png'); %read rgb image (wikipedia)
imtool(RGB);

BW = imcomplement(im2bw(RGB, graythresh(RGB))); %rgb -> bw + invert
% remove big white areas (percentage of white pixels in rectangular space
% of object high)
CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
imagesize = CC.ImageSize(1) * CC.ImageSize(2);
for i = 1:length(numPixels)
    x = floor(CC.PixelIdxList{i}/CC.ImageSize(1));
    y = mod(CC.PixelIdxList{i},CC.ImageSize(1));
    if numPixels(i)/((max(x)-min(x))*(max(y)-min(y))) > 0.167 && ((max(x)-min(x))*(max(y)-min(y)))/imagesize > 0.1
        BW(CC.PixelIdxList{i}) = 0;
    end
end
 imtool(BW);

BWT=bwmorph(BW, 'thin', inf);
imtool(BWT);

[H, T, R] = hough(BWT);

imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,...
      'InitialMagnification','fit');
title('Hough Transform');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(hot);
NHS=floor(((size(H)/25)-1)/2)*2+1; %NHoodSize Parameter
P  = houghpeaks(H,99,'threshold',ceil(0.1*max(H(:))),'NHoodSize',NHS);
x = T(P(:,2)); y = R(P(:,1));
plot(x,y,'s','color','white');
% Find lines and plot them
FG = length(BW); %FillGap Parameter
lines = houghlines(BWT,T,R,P,'FillGap',FG,'MinLength',7);
figure, imshow(RGB), hold on
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

end
