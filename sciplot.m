function [ ] = sciplot(Xdata, Ydata, varargin)
%% 本函数根据用户输入的数据实现快速绘图
%****************************************************************
%----             Author(s): Jia Yingqi                      ----
%----             Affiliation: Tongji University             ----
%----             E-mail: jiayingqi@tongji.edu.cn            ----
%----             Date: 06/02/2020                           ----
%****************************************************************
% 应用示例
% x1 = -0.5*pi:0.1:6.5*pi;
% y1 = 0.67*(1.4*cos(x1)+0.45)+5*rand(1,length(x1));
% x2 = -pi:0.01:5*pi;
% y2 = 1.45*sin(x2)+3.77;
% 
% green = [56 194 93]/256;
% blue = [76 114 176]/256;
% color = {green,blue};
% 
% sciplot({x1 x2}, {y1 y2}, 'color', color, 'linewidth', 3, 'fontsize', 25, ...
%     'xlabel', 'Time \itt \rm(s)', 'ylabel', 'Disp. \itx \rm(m)', ...
%     'sizeFigure', [100 100 800 500], 'legend', {'Eurocode values','Fitting values'}, ...
%     'grid', 'on', 'legendBox', 'off', 'save', './effect of sciplot')

%% 是否定义图窗大小
if isempty(find(strcmp(varargin,'sizeFigure'))) % 没有定义图窗大小
    figure('color', [1 1 1])
else % 定义了图窗大小
    sizeFigureValue = varargin{find(strcmp(varargin,'sizeFigure'))+1};
    figure('color', [1 1 1], 'position', sizeFigureValue)
end
    
%% 是否定义线条颜色
if isempty(find(strcmp(varargin,'color'))) % 如果没定义color，采用默认颜色
    blue = [76 114 176]/256;
    orange = [255 160 65]/256;
    purple = [184 135 195]/256;
    green = [56 194 93]/256;
    pink = [255 91 78]/256;
    gray = [164 160 155]/256;
    brown = [140 86 75]/256;
    yellow = [194 189 44]/256;
    colorValue = {blue orange purple green pink gray brown yellow};
else % 如果定义了color，采用定义的颜色
    colorValue = varargin{find(strcmp(varargin,'color'))+1};
end

%% 是否定义线条宽度
if isempty(find(strcmp(varargin,'linewidth'))) % 如果没定义color，采用默认线宽
    lineWidthValue = 2;
else 
    lineWidthValue = varargin{find(strcmp(varargin,'linewidth'))+1};
end

%% 绘制曲线
numCurve = length(Xdata); % 曲线的数量
for iCurve = 1:numCurve
    hold on
    plot(Xdata{iCurve}, Ydata{iCurve}, 'color', colorValue{iCurve}, ...
        'linewidth', lineWidthValue) 
end

%% 是否定义字号与字体
if isempty(find(strcmp(varargin,'fontname'))) % 如果没定义linewidth，采用默认字号
    fontNameValue = 'Times New Roman';
else
    fontNameValue = varargin{find(strcmp(varargin,'fontname'))+1};
end
if isempty(find(strcmp(varargin,'fontsize'))) % 如果没定义linewidth，采用默认字号
    fontSizeValue = 15;
else
    fontSizeValue = varargin{find(strcmp(varargin,'fontsize'))+1};
end
set(gca, 'fontname', fontNameValue, 'fontsize', fontSizeValue)

%% 是否定义标题与图例
if ~isempty(find(strcmp(varargin,'xlabel'))) % 如果定义x轴标题
    xLabelValue = varargin{find(strcmp(varargin,'xlabel'))+1};
    set(xlabel(xLabelValue), 'fontname', fontNameValue, 'fontsize', fontSizeValue)
end
if ~isempty(find(strcmp(varargin,'ylabel'))) % 如果定义y轴标题
    yLabelValue = varargin{find(strcmp(varargin,'ylabel'))+1};
    set(ylabel(yLabelValue), 'fontname', fontNameValue, 'fontsize', fontSizeValue)
end
if ~isempty(find(strcmp(varargin,'title'))) % 如果定义图像标题
    titleValue = varargin{find(strcmp(varargin,'title'))+1};
    set(title(titleValue), 'fontname', fontNameValue, 'fontsize', fontSizeValue)
end
if ~isempty(find(strcmp(varargin,'legend'))) % 如果定义图例
    legendValue = varargin{find(strcmp(varargin,'legend'))+1};
    gray = [164 160 155]/256;
    set(legend(legendValue), 'fontname', fontNameValue, 'fontsize', fontSizeValue, ...
        'EdgeColor', gray, 'linewidth', 1.5)
end

%% 坐标显示定义
% x轴坐标显示
lb = min(min(cell2mat(Xdata)));
ub = max(max(cell2mat(Xdata)));
[axis_display_data, int_digit, decimal_digit] = axis_boundary(lb, ub);
xlim([axis_display_data(1) axis_display_data(end)]) % 坐标范围修改
set(gca, 'xtick', axis_display_data) % 坐标显示值修改
set(gca, 'xticklabel', sprintf(['%.',num2str(decimal_digit),'f\n'], ...
    get(gca,'xtick'))) % 坐标精度修改
% y轴坐标显示
lb = min(min(cell2mat(Ydata)));
ub = max(max(cell2mat(Ydata)));
[axis_display_data, int_digit, decimal_digit] = axis_boundary(lb, ub);
ylim([axis_display_data(1) axis_display_data(end)]) % 坐标范围修改
set(gca, 'ytick', axis_display_data) % 坐标显示值修改
set(gca, 'yticklabel', sprintf(['%.',num2str(decimal_digit),'f\n'], ...
    get(gca,'ytick'))) % 坐标精度修改

%% 其他定义
box on
set(gca, 'looseInset', [0 0 0 0], 'linewidth', 2)

%% 是否开启网格（默认开启）
set(gca,'GridLineStyle',':','GridColor','k')
grid on
if ~isempty(find(strcmp(varargin,'grid')))
    flag = varargin{find(strcmp(varargin,'grid'))+1};
    if strcmp(flag, 'off') 
        grid off
    end
end

%% 文件保存定义
if ~isempty(find(strcmp(varargin,'save')))
    fileName = varargin{find(strcmp(varargin,'save'))+1};
    print(fileName, '-djpeg', '-r300');
end

function [axis_display_data, int_digit, decimal_digit] = axis_boundary(lb,ub)
%% 本函数根据输入数据确定合理的坐标显示范围
% 应该根据data的极差，将其分为4-6份，每份的大小（位数）应该为最大数据的位数（或位数减一）
% lb = min(data) 数据下界
% ub = max(data) 数据上界

%% 分析数据
if lb*ub <= 0 % 如果跨越坐标轴原点
    flag = 1; % 根据绘图规范，坐标显示范围必须包含原点
else
    flag = 0;
end
range = ub-lb; % 数据的极差
maxData = max(abs([lb,ub])); % 绝对值最大的数据

%% 将极差化为整数
stringData = num2str(maxData); % 将最大数据maxData化为字符串
if ~isempty(find(stringData=='.')) % 判断maxData是否为小数，小数必然存在小数点
    decimal_digit = length(stringData) - find(stringData=='.'); % 小数部分的位数
    if stringData(1) == '0' % 纯小数（pure decimal）
        all_digit = length(stringData)-2; % 化为整数后的全部位数
    else % 带小数（mixed decimal）
        all_digit = length(stringData)-1; % 化为整数后的全部位数
    end
    amp = decimal_digit;
    temp_range = 10^amp*range; % 将极差化为整数
    temp_digit = floor(log10(temp_range))+1; % 新极差的位数
    % 对于任意正实数N，若其为m位数，则有：10^(m-1)<=N<10^m
    % 从而有:m-1<=log10(N)<m → log10(N)<m<=log10(N)+1且m为整数 
    % 因此，m为log10(N)+1的向下取整
    if temp_digit == 1
        amp = amp+1; % 对于个位数，需要化为十位数
    end
else % maxData为整数
    digit = floor(log10(range))+1; % 极差的位数
    all_digit = length(num2str(maxData)); % 化为整数后的全部位数
    if digit == 1
        amp = 1; % 对于个位数，需要化为十位数
    else
        amp = 0;
    end
end
temp_range = round(10^amp*range); % 由于精度问题，必须使用round
% 每份的位数只能为m或m-1，
% 即-20~20可以被分为-20:5:20，也可以被分为-20:10:20
% 因此，需要对其暂时乘以10^amp后再进行如下判断

%% 每份的大小需要被divisor_delta整除
if all_digit >= 3
    amp1 = all_digit-2;
    divisor_delta = 5*10^amp1; 
else
    divisor_delta = 1;
end

%% 确定坐标显示范围
axis_display_data = [];
for divisor = [4,5,6] % 允许将数据分为的份数
    if rem(temp_range,divisor) == 0 % 数据能够等分为若干份
        temp_axis_display_data = linspace(lb,ub,divisor+1); % 拟定的坐标显示范围
        delta = temp_range/divisor; % 每份的大小
        if rem(delta,divisor_delta) == 0 % 每份的大小能够被divisor_delta整除
            delta = delta/10^amp; % delta的真实值
            if flag == 1 % 坐标显示范围有必要包含原点
                if ~isempty(find(temp_axis_display_data==0)) % 当前拟定范围包含原点
                    axis_display_data = temp_axis_display_data;
                    break
                end
            else % 没有必要包含原点 
                axis_display_data = temp_axis_display_data;
                break
            end
        end
    end
end

if isempty(axis_display_data) % 判断是否仍然为空集 
    divisor = 4; % 若仍然为空集，拟定分为4份
    delta = fix(temp_range/divisor);
    temp = round(delta/divisor_delta)*divisor_delta; % 每份的大小需要被divisor_delta整除
    if temp ~= 0
        delta = round(delta/divisor_delta)*divisor_delta;
    else
        divisor_delta = divisor_delta/5;
        delta = round(delta/divisor_delta)*divisor_delta;
    end
    delta = delta/10^amp; % delta的真实值
    if flag == 1 % 坐标显示范围有必要包含原点
        num_ub = ceil(ub/delta); % 原点以上的份数
        num_lb = ceil(-lb/delta); % 原点以下的份数
        axis_display_data = -num_lb*delta:delta:num_ub*delta;
    else % 没有必要包含原点
        num = ceil((ub-lb)/delta); % 份数
        lb = floor(10^amp*lb/divisor_delta)*divisor_delta/10^amp;
        if lb+num*delta >= ub
            axis_display_data = linspace(lb, lb+num*delta, num+1);
        else
            axis_display_data = linspace(lb, lb+(num+1)*delta, num+2);
        end
    end
end

%% 确定坐标显示精度
stringData = num2str(delta); % 将每份的大小delta化为字符串
if ~isempty(find(stringData=='.')) % 判断delta是否为小数，小数必然存在小数点
    decimal_digit = length(stringData) - find(stringData=='.'); % 小数部分的位数
    int_digit = find(stringData=='.')-1; % 整数部分的位数
else
    decimal_digit = 0;
    int_digit = length(stringData);
end