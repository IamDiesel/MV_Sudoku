% mvNextSlide

% Automatische Berechnung der Foliennummer
fn = fn + 1; fntext = sprintf('Folie Nr. %d',fn);

% figure- und plot-Fenster definierter Größe 

axes('Position', plotKoord,'Visible','on','box','on'); 
set(f, 'MenuBar', 'none');
set(f, 'ToolBar', 'none');


% Ausgabe der Kopfzeile
text(0.015, kopfzeilpos, mvKopfzeile,'FontUnits','normalized','FontSize',fkopf);

% Ausgabe des Folientitels
text(linkerEinzug, titelpos, mvTitel,'Interpreter','latex','FontUnits','normalized','FontSize',ftitel);

% Ausgabe des Folientextes
text(linkerEinzug, textpos, mvText ,'Interpreter','latex','FontUnits','normalized','FontSize',ftext); 

% Ausgabe der drei Fusszeilentexte
text(0.015, fusszeilpos, mvFusszeileL,'FontUnits','normalized','FontSize',ffuss); 
text(0.47,  fusszeilpos, fntext,'FontUnits','normalized','FontSize',ffuss);
text(0.75,  fusszeilpos, nextslide,'FontUnits','normalized','FontSize',ffuss);