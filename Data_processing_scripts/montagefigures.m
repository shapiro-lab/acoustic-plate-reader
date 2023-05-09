function layout = montagefigures(figs_array,n_rows,n_cols)
% figurearray is a gobjects array for the figures that you'd like
% included in the montage
layout = tiledlayout(n_rows,n_cols)
for wellIx = 1:length(figs_array)
    curr_axes = gca;
    disp(wellIx);
    fig = figs_array(wellIx).CurrentAxes
    graphics = fig.Children; %CurrentObject
    for graph_index = 1:length(graphics)
        copyobj(graphics(graph_index),curr_axes,'legacy');
    end
    nexttile
end
end
