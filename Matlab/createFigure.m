function f = createFigure(figureSize)
c = [...
    0.6   1      0.7    ;...
    1      1  0.9   ;...
    0.1  1 1      ;...
    0.8 1 0.6    ;...
    0.3333  1 0.5;...
    0.58   0.67   0.9];

screenSize = get(0,'ScreenSize');

f = figure;

% set line properties:
set(f,'DefaultLineLineWidth',1.5);
set(f,'DefaultLineMarkerSize',4);

% set axes Properties
set(f,'DefaultAxesLineWidth',1);
%set(f,'DefaultAxesFontSize',12);
set(f,'DefaultAxesFontName','Arial');
set(f,'DefaultAxesColor','None');
set(f,'DefaultAxesColorOrder',hsv2rgb(c));

% set figure properties
set(f,'Color','White');
figurePos = get(f,'Position');
set(f,'Position',...
    [figurePos(1) screenSize(4)-figurePos(2),...
    figurePos(3:4)]);
end