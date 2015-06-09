
imtool close all;
close all;
clear;
mvImageClose = 0.1;
mvInitSlides % mvTitel und mvText sind noch leer ...

%% -----------------------   Folie    ------------------------------------
mvTitel = {
    ''
    ''
    ''
    '\bf{MV Projekt: Sudoku-Zahlenerkennung}'
    };

mvText = {
    'Optische Zahlenerkennung in Sudokufeldern mittels Matlab'
    ''
    ''
    ''
    'von'
    ''
    mvAutor1
    mvAutor2
    mvAutor3
    ''
    ''
    ''
    ''
    };
mvOpenSlide;  % Die vorbereiteten Texte werden geschrieben.
mvClearSlide; %mvCloseSlide % pause(); Tastendruck schliesst Folie und löscht mvTitel und mvText.
% ------------------------------------------------------------------------

%% -----------------------   Folie    ------------------------------------
mvTitel = '\bf{Prozesskette Soduko-Zahlenerkennung}';
mvText = {
    'Laden des Bildes'
    ''
    'Grauwertkonvertierung und Schwellwertbestimmung'
    ''
    'Binärkonvertierung & Invertierung'
    ''
    'Entfernung von Schatten'
    ''
    'Sudokufeld-Randerkennung'
    ''
    'Sudokufeld-Transformation'
    ''
    'Erkennen der Zahlenfelder'
    ''
    'Merkmalsextraktion'
    ''
    'Zahlenerkennung (Neuronales Netz)'
     };
mvNextSlide%mvOpenSlide  % Die vorbereiteten Texte werden geschrieben.
mvClearSlide; %mvCloseSlide % pause(); Tastendruck schliesst Folie und löscht mvTitel und mvText.
% ------------------------------------------------------------------------

%% -----------------------   Folie Vorbereitungen  ------------------------------------
mvNextSlide%mvOpenSlide; % Leere Folienfläche
text(0.01, 0.9,'-Laden des Bildes mittels:','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
text(0.1, 0.82, 'I = imread(''sudoku6.jpg'');','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0.1 0.1 0.1]);
pause();
        I = imread('sudoku6.jpg');
        hFig = imtool(I); %show orignal Image
        pause(mvImageClose);
        close(hFig);      %close original Image after 10s
text(0.01, 0.7, '-Grauwertkonvertierung mittels rgb2gray','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
        if numel(size(I)) == 3
            IG = rgb2gray(I);
        else
            IG = I;
        end
        imageArea = numel(IG); % number of pixels
        imageSize = size(IG); % width and height
pause();
text(0.01, 0.6, '-Bestimmung des Schwellwertes für die Binärbildumwandlung:','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
text(0.1, 0.52, 'grayThreshold = median(double(IG(:)))/1.2/255;','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0.1 0.1 0.1]);
        grayThreshold = median(double(IG(:)))/1.2/255; % grayThreshold
pause();
text(0.01, 0.4, '-Umwandlung in Bin\"arbild und Invertierung','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
text(0.1, 0.32, 'BW = imcomplement(im2bw(IG,grayThreshold));','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0.1 0.1 0.1]);
        pause();
        BW = imcomplement(im2bw(IG,grayThreshold)); %gray -> bw + invert
        hFig = imtool(BW);% black and white image
        pause(mvImageClose);
        close(hFig);       
        pause();
 text(0.01, 0.2, '-Entfernung von Schatten am Rand: ','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
 text(0.1, 0.12, 'BW = imclearborder(BW);','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0.1 0.1 0.1]);
        pause();
        BW = imclearborder(BW);
        hFig = imtool(BW);% black and white image
        pause(mvImageClose);
        close(hFig);       
        pause();
mvClearSlide; %mvCloseSlide % pause(); Tastendruck schliesst Folie und löscht mvTitel und mvText.
% ------------------------------------------------------------------------

%% -----------------------   Folie Sudokufeld Randerkennung    ------------------------------------
mvNextSlide%mvOpenSlide  % Die vorbereiteten Texte werden geschrieben.
text(0.01, 0.9,'-Regionprops erhalten:','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
text(0.1, 0.82, 'R = regionprops(BW,''Area'',''BoundingBox'',''PixelList'');','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0.1 0.1 0.1]);
            pause();
text(0.01, 0.72,'-Iteration über die Regionprops R','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
text(0.01, 0.62,'-->Größtes Regionprop finden','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
            pause();
text(0.01, 0.52,'-Ecken des größten Regionprobs berechnen','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
            pause();
            
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
mvClearSlide; %mvCloseSlide % pause(); Tastendruck schliesst Folie und löscht mvTitel und mvText.
% ------------------------------------------------------------------------

%% -----------------------   Folie Sudokufeld Transformation  ------------------------------------
mvNextSlide%mvOpenSlide; % Leere Folienfläche
text(0.01, 0.9,'-Erstellen der Transformation und Anwendung der Transformation auf das Bild:','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
text(0.1, 0.82, 'T = fitgeotrans(pts(1:4,:),900*[0 0; 1 0; 1 1; 0 1],''projective'');','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0.1 0.1 0.1]);
text(0.1, 0.74, 'R = regionprops(IT,''Area'',''BoundingBox'');','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0.1 0.1 0.1]);
pause();
        T = fitgeotrans(pts(1:4,:),900*[0 0; 1 0; 1 1; 0 1],'projective');
        IT = logical(imwarp(double(BW),T));
        hFig = imtool(IT); %show orignal Image
        pause(mvImageClose);
        close(hFig);      %close original Image after 10s
pause();
text(0.01, 0.64,'-In Transformiertem Bild größtes Regionprop finden (vgl. Iteration Randerkennung)','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
        R = regionprops(IT,'Area','BoundingBox');
        NR = numel(R);
        kmax = 0;
        maxArea = 0;
        for k = 1:NR
            if R(k).Area > maxArea
                maxArea = R(k).Area;
                kmax = k;
            end
        end
        BBmax = R(kmax).BoundingBox;
pause();
text(0.01, 0.54,'-Sudoku ausschneiden und als eigenständiges Bild betrachten','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
text(0.1, 0.46,'Hilfslinien verdeutlichen die gefundenen Elemente des Bildes','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
        IT = IT(ceil(BBmax(2)):ceil(BBmax(2)+BBmax(4)),ceil(BBmax(1)):ceil(BBmax(1)+BBmax(3)));
        IT = bwmorph(IT,'close');
        imageArea = numel(IT); % number of pixels
        imageSize = size(IT); % width and height

        %% remove noise
        IT = bwareaopen(IT,100); % opening

        R = regionprops(IT, 'BoundingBox', 'Area', 'PixelIdxList');
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
        
        IRahmen = imread('sudoku6_Rahmen.tif');
        hFig = imtool(IRahmen);
        pause(mvImageClose);
        close(hFig);
        pause();


mvClearSlide; %mvCloseSlide % pause(); Tastendruck schliesst Folie und löscht mvTitel und mvText.
% ------------------------------------------------------------------------

%% -----------------------   Folie Erkennen der Zahlenfelder   ------------------------------------
mvNextSlide%mvOpenSlide; % Leere Folienfläche
text(0.01, 0.9,'-Zahlenfelder erkennen mit Hilfe von Area-Eigenschaften:','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
text(0.1, 0.82, 'Vergleich der Area-Eigenschaft aller Elemente im Bild','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0.1 0.1 0.1]);
text(0.1, 0.74, 'Felder sind größer wie Pixel/3000 und kleiner Pixel/90','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0.1 0.1 0.1]);
pause();
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
            
            
            Pnew(k,:) = [R(k).BoundingBox(1)+R(k).BoundingBox(3)/2 R(k).BoundingBox(2)+R(k).BoundingBox(4)/2];       
            if inpolygon(Pnew(k,1),Pnew(k,2),pts(:,1),pts(:,2))
                %h_digitcircles(k) = plot(Pnew(k,1),Pnew(k,2),'ro','markersize',max([R(k).BoundingBox(3:4) 20]));
                kgood(k) = 1;
            end      
    end
end

text(0.01, 0.66, '-Merkmalsextraktion durch Auslesen der Regionprops aller gefundenen Elemente','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0.1 0.1 0.1]);
        %% generate input
        R = regionprops(IT,'all');
        R = R(2:size(R));

        for i = 1:size(R)
                    input(1,i) = (R(i).Centroid(1)-R(i).BoundingBox(1))/R(i).BoundingBox(3);
                    input(2,i) = (R(i).Centroid(2)-R(i).BoundingBox(2))/R(i).BoundingBox(4);
                    input(3,i) = R(i).MajorAxisLength/sqrt(sum(R(i).BoundingBox(3:4).^2));
                    input(4,i) = R(i).MinorAxisLength/sqrt(sum(R(i).BoundingBox(3:4).^2));
                    input(5,i) = R(i).Eccentricity;
                    input(6,i) = R(i).Orientation;
                    input(7,i) = R(i).ConvexArea/(prod(R(i).BoundingBox(3:4)));
                    input(8,i) = R(i).FilledArea/(prod(R(i).BoundingBox(3:4)));
                    input(9,i) = R(i).EulerNumber;
                    input(10,i) = R(i).Solidity;
                    input(11,i) = R(i).Extent;
                    input(12,i) = R(i).Perimeter;

        end
text(0.01, 0.58, '-Übergabe der Merkmale in bereits gelerntes Neuronales Netz','Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0.1 0.1 0.1]);
        load('net.mat');
        numbers = net(input);
        [temp, blub] = max(numbers);
        numbersI = round(temp).*blub;
        sudoku = zeros(9,9);
        for i=1:size(R)
            sudoku(floor(R(i).BoundingBox(2)/100)+1,floor(R(i).BoundingBox(1)/100)+1) = numbersI(i);
        end
        
pause();
mvClearSlide; %mvCloseSlide % pause(); Tastendruck schliesst Folie und löscht mvTitel und mvText.
% ------------------------------------------------------------------------

% %% -----------------------   Folie Ausgabe gefundene Zahlen   ------------------------------------
mvNextSlide%mvOpenSlide; % Leere Folienfläche
  Sudoku_Found = cell(10, 10);
%Erstellen der Kopfzeile und der 1. Spalte
for i = 1:9
    Sudoku_Found{1, i+1} = sprintf('Col%d', i);     % Not recommended
    Sudoku_Found{i+1,1} = sprintf('Row%d', i);
end
%Übertragen der gefundenen Werte in das Cell-Array
for i=2:10
    for j =2:10
    Sudoku_Found{i,j} = sudoku(i-1,j-1);
    end
end

%Ausgabe des Cell-Aray in Folie
for i=2:10
    dec = i-1;
   text(0.1*dec, 0.95,Sudoku_Found{1,i},'Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]); 
   text(0.01, 1-0.1*dec,Sudoku_Found{i,1},'Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]); 
end

 
for i=2:10
    dec_i = i-1;
    for j=2:10
        dec_j = j-1;
        if Sudoku_Found{j,i} ~= 0
            text(0.1*dec_i, 1-0.1*dec_j,num2str(Sudoku_Found{j,i}),'Interpreter','latex','FontUnits','normalized','FontSize',0.05, 'Rotation',0,'Color',[0 0 0]);
        end
    end
end

pause();
mvClearSlide; %mvCloseSlide % pause(); Tastendruck schliesst Folie und löscht mvTitel und mvText.
% ------------------------------------------------------------------------

%% -----------------------   Folie    ------------------------------------
mvTitel = '\bf{Vielen Dank f\"ur Ihre Aufmerksamkeit.}';        
mvNextSlide%mvOpenSlide;
mvCloseSlide % pause(); Tastendruck schliesst Folie und löscht mvTitel und mvText.
% ------------------------------------------------------------------------