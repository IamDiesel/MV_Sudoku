imtool close all;
close all;
clear;
RGB = imread('../sudoku8.jpg'); %read rgb image (photo)
%RGB = imread('../sudoku1.png'); %read rgb image (wikipedia)
imtool(RGB);

BW = imcomplement(im2bw(RGB,graythresh(RGB))); %rgb -> bw + invert
imtool(BW);
% remove big white areas (percentage of white pixels in rectangular space
% of object high)
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

BWT=bwmorph(BW, 'thin', inf);
imtool(BWT);

[H, T, R] = hough(BWT);

imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,...
      'InitialMagnification','fit');
title('Hough Transform');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(hot);
NHS=floor(((size(H)/20)-1)/2)*2+1; %NHoodSize Parameter
P  = houghpeaks(H,99,'threshold',ceil(0.3*max(H(:))),'NHoodSize',NHS);
x = T(P(:,2)); y = R(P(:,1));
plot(x,y,'s','color','white');
% Find lines and plot them
FG = length(BW); %FillGap Parameter
lines = houghlines(BWT,T,R,P,'FillGap',FG,'MinLength',7);
% figure, imshow(RGB), hold on
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
%    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
%    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

end

figure, imshow(RGB), hold on
for i = 1:length(lines)
   line(i,1) = lines(i).point1(1);
   line(i,2) = lines(i).point1(2);
   line(i,3) = lines(i).point2(1) - lines(i).point1(1);
   line(i,4) = lines(i).point2(2) - lines(i).point1(2);
   plot([line(i,1) line(i,1)+line(i,3)],[line(i,2) line(i,2)+line(i,4)],'LineWidth',2,'Color','green')
end
n=1;
for i = 1:length(lines)
    for j = i+1:length(lines)
       r1  = line(i,3)/line(j,3);
       r2 = line(i,4)/line(j,4);
       if (r1 ~= r2 && i ~=j)
           A = [line(i,3) -line(j,3);
                line(i,4) -line(j,4)]; 
           b = [line(j,1)-line(i,1);
               line(j,2)-line(i,2)];
           temp = A\b;
           if ((abs(temp(1)) ~= Inf) || (abs(temp(1)) ~= Inf))
                L(:,n) = [i;
                          j;
                          line(i,1)+temp(1)*line(i,3);
                          line(i,2)+temp(1)*line(i,4)];
                plot(L(3,n),L(4,n),'x','LineWidth',2,'Color','red');
                n = n+1;
           end
       end
    end
end

