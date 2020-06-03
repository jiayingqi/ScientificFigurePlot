clear
clc
close all

%% sciplot”¶”√ æ¿˝
x1 = -0.5*pi:0.1:6.5*pi;
y1 = 0.67*(1.4*cos(x1)+0.45)+5*rand(1,length(x1));
x2 = -pi:0.01:5*pi;
y2 = 1.45*sin(x2)+3.77;

green = [56 194 93]/256;
blue = [76 114 176]/256;
color = {green,blue};

sciplot({x1 x2}, {y1 y2}, 'color', color, 'linewidth', 3, 'fontsize', 25, ...
    'xlabel', 'Time \itt \rm(s)', 'ylabel', 'Disp. \itx \rm(m)', ...
    'sizeFigure', [100 100 800 500], 'legend', {'Eurocode values','Fitting values'}, ...
    'grid', 'on', 'legendBox', 'off', 'save', './effect of sciplot')