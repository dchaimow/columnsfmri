function sim_plotPattern(sim,y,titleStr,cmap)
if ~exist('cmap','var')
    cmap = 'gray';
end
if ~exist('titleStr','var')
    titleStr = '';
end
minx = min(sim.x)+sim.L/(2*size(y,1));
maxx = max(sim.x)-sim.L/(2*size(y,2));

imagesc([minx maxx],[minx maxx],y);
title(titleStr);
colormap(gca,cmap);
colorbar;
xlabel('position [mm]');
axis image;
set(gca,'FontSize',14);

