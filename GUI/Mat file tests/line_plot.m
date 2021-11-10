fig = figure;
hold on
for k = 1:10
%   When you create a plot, you can specify the legend labels by setting the “DisplayName” property as name-value pair.
%   Set the "DisplayName" property to a character vector of the text that you want to include in the legend.
%   To include a variable value in the text, use “num2str”.
    txt = ['X = ',num2str(k)];
    plot(rand(10,1), rand(10,1),'DisplayName',txt)
end
hold off
legend show