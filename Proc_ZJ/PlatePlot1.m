function fig_out = PlatePlot1(data,fig_title,fontsize,varargin)
if nargin == 2
    fontsize = 16;
elseif nargin == 1
    fig_title = '';
    fontsize = 16;
end
fig_out = figure;
[y,x] = size(data);
imagesc(data);axis image
xticks(1:x)
yticks(1:y)
ylabelnames = {'A','B','C','D','E','F','G','H'};
ylabelnames = ylabelnames(1:y);
yticklabels(ylabelnames)
colorbar;
title(fig_title)
set(gca,'fontsize',fontsize);