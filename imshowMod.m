subplot(2,3,1); imshow(image); colormap 'gray';
title(filename, 'FontUnits','normalized','FontSize',ftitel);
freezeColors;

subplot(2,3,2); imshow(image);%colormap 'jet';
colorbar();
title(filename, 'FontUnits','normalized','FontSize',ftitel);

% Histogramm des Graubildes anzeigen
subplot(2,3,3); imhist(image);
ti = title('Histogramm', 'FontUnits','normalized','FontSize',ftitel);

[srad, sang, S] = specxture(image);
% Spektralbild anzeigen
subplot(2,3,4); imshow(S); colormap 'jet';
title('Ortsfrequenzspektrum S', 'FontUnits','normalized','FontSize',ftitel);
%
subplot(2,3,5); plot(srad);
title('Radiale Spektralfunktion srad', 'FontUnits','normalized','FontSize',ftitel);
%
subplot(2,3,6); plot(sang);
title('Angulare Spektralfunktion sang', 'FontUnits','normalized','FontSize',ftitel);

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.05,'\bf Abbrechen mit ctrl-C, Weiter mit Leertaste','Interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'top','FontUnits','normalized','FontSize',.02);
text(0.5, 0.99 ,mvtitel,'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment', 'top','FontUnits','normalized','FontSize',.03);
