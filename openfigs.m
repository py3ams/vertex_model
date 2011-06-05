function openfigs(fig_saves_name)

fig_saves_location = ['Figs/',fig_saves_name,'/'];

if ~exist(fig_saves_location,'dir')
    error([fig_saves_location,' does not exist'])
end

h = dir([fig_saves_location,'*.fig']);

no_figs = length(h);

for i = 1:no_figs
    open([fig_saves_location,h(i).name]) 
end