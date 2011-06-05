disp('busy');

disp('Click on subplot you wish to copy and hit return')
pause;

ax1 = gca;

% create a new figure
hf = figure('position',[100 100 300 300],'PaperPositionMode','auto');
% create a dummy axes in the standard position
dax = axes;
% store the position of the dummy axes
pos = get(dax,'position');
pos = [pos(1)+0.03 pos(2)+0.02 pos(3) pos(4)];
% delete the dummy axes
delete(dax);

% copy the subplot into the new figure
ax2 = copyobj(ax1,hf);
% resize the axes to take up the whole figure
set(ax2,'position',pos);

% simulation_name = 'faster_lateral_growth';

% set(gca,'XTickLabel',sprintf('%0.0f|',str2num(get(gca,'XTickLabel'))))

% set(gca,'YTick',[0:0.1:1]*10^-4)
% set(gca,'YTickLabel',sprintf('%0.2f|',str2num(get(gca,'YTickLabel'))))

% axis([0 10*10^4 0 1.1*10^-4])

title('')
xlabel('Time','FontName','arial')
ylabel('Chemical','FontName','arial')

set(gca,'FontName','arial');

% saveas(gcf,['Figs/',simulation_name,'/mitosis_locations.eps'],'psc2')




