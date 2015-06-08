% mvInitSlides
% Vorlage f�r MATLAB-Pr�sentationsfolien
% Reinhard Malz, SS 2015
% 2015-04-28 Version 2.0

%% Eingabebereich f�r alle Folien

% Kopfzeile mit projektspezifischen Angaben 
mvProjekt = '"Sudoku Zahlenerkennung"';
mvAutor1  = 'Daniel Kahrizi';
mvAutor2  = 'Nikolas Porzer';
mvAutor3  = 'Matthias Roth';
mvAutor4  = '';
mvKopfzeile = sprintf('Projekt:  %s  von  %s  %s  %s   %s', ...
              mvProjekt, mvAutor1, mvAutor2, mvAutor3, mvAutor4);

% Folientitel und Folientext 
mvTitel = '';
mvText = '';
% Bleiben diese Variablen leer, �ffnet mvOpenSlide eine leere Folie, die 
% mit eigenen Formatierungen beschrieben werden kann. Siehe mvFolienDemo.m 

% Fusszeile  
nextslide = '>>> Weiter mit Leertaste oder Esc >>>';
mvFusszeileL = sprintf('Maschinelles Sehen SS 2015');

% normalisierte Fontgr��en
ftitel = 0.04;% 0.05;
ftext = 0.03;%0.04;
ffuss = 0.025;
fkopf = 0.025;

% Vollformat f�r graues figure Fenster
figuKoord = [0.02 0.06 0.96 0.88];

if 0
    % Wei�es Plotfenster vollformatig
    plotKoord = [0 0 .95 .95];
    % Standardformatierungen bezogen auf wei�es Plotfenster
    kopfzeilpos  = 0.975; % 1.0 ist oberer Fensterrand
    titelpos = 0.85;
    textpos = 0.4;
    fusszeilpos =  0.025; % 0.0 ist unterer Fensterrand
else
    % Wei�es Plotfenster mit Platz f�r Kopf- und Fusszeile
    plotKoord = [0.02 0.06 0.96 0.88];
    % Standardformatierungen bezogen auf wei�es Plotfenster
    kopfzeilpos  = 1.025; % 1.0 ist oberer Fensterrand
    titelpos = 0.85;
    textpos = 0.4;
    fusszeilpos = -0.04; % 0.0 ist unterer Fensterrand
end

linkerEinzug = 0.05; % etwas einger�ckt

fn = 0;  % Foliennummer
