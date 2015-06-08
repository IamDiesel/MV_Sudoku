
imtool close all;
close all;
clear;
mvImageClose = 2;
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
mvClearSlide; %mvCloseSlide % pause(); Tastendruck schliesst Folie und löscht mvTitel und mvText.
% ------------------------------------------------------------------------

%% -----------------------   Folie    ------------------------------------
mvNextSlide%mvOpenSlide; % Leere Folienfläche
imshow(imread('lena.tif'));
pause(1);
text(10,20, 'Ich bin kein Roboter!','Interpreter','latex','FontUnits','normalized','FontSize',0.1, 'Rotation',0,'Color',[1 1 0]);
mvClearSlide; %mvCloseSlide % pause(); Tastendruck schliesst Folie und löscht mvTitel und mvText.
% ------------------------------------------------------------------------

%% -----------------------   Folie    ------------------------------------
mvTitel = '\bf{Vielen Dank f\"ur Ihre Aufmerksamkeit.}';        
mvNextSlide%mvOpenSlide;
mvCloseSlide % pause(); Tastendruck schliesst Folie und löscht mvTitel und mvText.
% ------------------------------------------------------------------------