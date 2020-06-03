function [ ] = sciplot(Xdata, Ydata, varargin)
%% �����������û����������ʵ�ֿ��ٻ�ͼ
%****************************************************************
%----             Author(s): Jia Yingqi                      ----
%----             Affiliation: Tongji University             ----
%----             E-mail: jiayingqi@tongji.edu.cn            ----
%----             Date: 06/02/2020                           ----
%****************************************************************
% Ӧ��ʾ��
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

%% �Ƿ���ͼ����С
if isempty(find(strcmp(varargin,'sizeFigure'))) % û�ж���ͼ����С
    figure('color', [1 1 1])
else % ������ͼ����С
    sizeFigureValue = varargin{find(strcmp(varargin,'sizeFigure'))+1};
    figure('color', [1 1 1], 'position', sizeFigureValue)
end
    
%% �Ƿ���������ɫ
if isempty(find(strcmp(varargin,'color'))) % ���û����color������Ĭ����ɫ
    blue = [76 114 176]/256;
    orange = [255 160 65]/256;
    purple = [184 135 195]/256;
    green = [56 194 93]/256;
    pink = [255 91 78]/256;
    gray = [164 160 155]/256;
    brown = [140 86 75]/256;
    yellow = [194 189 44]/256;
    colorValue = {blue orange purple green pink gray brown yellow};
else % ���������color�����ö������ɫ
    colorValue = varargin{find(strcmp(varargin,'color'))+1};
end

%% �Ƿ����������
if isempty(find(strcmp(varargin,'linewidth'))) % ���û����color������Ĭ���߿�
    lineWidthValue = 2;
else 
    lineWidthValue = varargin{find(strcmp(varargin,'linewidth'))+1};
end

%% ��������
numCurve = length(Xdata); % ���ߵ�����
for iCurve = 1:numCurve
    hold on
    plot(Xdata{iCurve}, Ydata{iCurve}, 'color', colorValue{iCurve}, ...
        'linewidth', lineWidthValue) 
end

%% �Ƿ����ֺ�������
if isempty(find(strcmp(varargin,'fontname'))) % ���û����linewidth������Ĭ���ֺ�
    fontNameValue = 'Times New Roman';
else
    fontNameValue = varargin{find(strcmp(varargin,'fontname'))+1};
end
if isempty(find(strcmp(varargin,'fontsize'))) % ���û����linewidth������Ĭ���ֺ�
    fontSizeValue = 15;
else
    fontSizeValue = varargin{find(strcmp(varargin,'fontsize'))+1};
end
set(gca, 'fontname', fontNameValue, 'fontsize', fontSizeValue)

%% �Ƿ��������ͼ��
if ~isempty(find(strcmp(varargin,'xlabel'))) % �������x�����
    xLabelValue = varargin{find(strcmp(varargin,'xlabel'))+1};
    set(xlabel(xLabelValue), 'fontname', fontNameValue, 'fontsize', fontSizeValue)
end
if ~isempty(find(strcmp(varargin,'ylabel'))) % �������y�����
    yLabelValue = varargin{find(strcmp(varargin,'ylabel'))+1};
    set(ylabel(yLabelValue), 'fontname', fontNameValue, 'fontsize', fontSizeValue)
end
if ~isempty(find(strcmp(varargin,'title'))) % �������ͼ�����
    titleValue = varargin{find(strcmp(varargin,'title'))+1};
    set(title(titleValue), 'fontname', fontNameValue, 'fontsize', fontSizeValue)
end
if ~isempty(find(strcmp(varargin,'legend'))) % �������ͼ��
    legendValue = varargin{find(strcmp(varargin,'legend'))+1};
    gray = [164 160 155]/256;
    set(legend(legendValue), 'fontname', fontNameValue, 'fontsize', fontSizeValue, ...
        'EdgeColor', gray, 'linewidth', 1.5)
end

%% ������ʾ����
% x��������ʾ
lb = min(min(cell2mat(Xdata)));
ub = max(max(cell2mat(Xdata)));
[axis_display_data, int_digit, decimal_digit] = axis_boundary(lb, ub);
xlim([axis_display_data(1) axis_display_data(end)]) % ���귶Χ�޸�
set(gca, 'xtick', axis_display_data) % ������ʾֵ�޸�
set(gca, 'xticklabel', sprintf(['%.',num2str(decimal_digit),'f\n'], ...
    get(gca,'xtick'))) % ���꾫���޸�
% y��������ʾ
lb = min(min(cell2mat(Ydata)));
ub = max(max(cell2mat(Ydata)));
[axis_display_data, int_digit, decimal_digit] = axis_boundary(lb, ub);
ylim([axis_display_data(1) axis_display_data(end)]) % ���귶Χ�޸�
set(gca, 'ytick', axis_display_data) % ������ʾֵ�޸�
set(gca, 'yticklabel', sprintf(['%.',num2str(decimal_digit),'f\n'], ...
    get(gca,'ytick'))) % ���꾫���޸�

%% ��������
box on
set(gca, 'looseInset', [0 0 0 0], 'linewidth', 2)

%% �Ƿ�������Ĭ�Ͽ�����
set(gca,'GridLineStyle',':','GridColor','k')
grid on
if ~isempty(find(strcmp(varargin,'grid')))
    flag = varargin{find(strcmp(varargin,'grid'))+1};
    if strcmp(flag, 'off') 
        grid off
    end
end

%% �ļ����涨��
if ~isempty(find(strcmp(varargin,'save')))
    fileName = varargin{find(strcmp(varargin,'save'))+1};
    print(fileName, '-djpeg', '-r300');
end

function [axis_display_data, int_digit, decimal_digit] = axis_boundary(lb,ub)
%% ������������������ȷ�������������ʾ��Χ
% Ӧ�ø���data�ļ�������Ϊ4-6�ݣ�ÿ�ݵĴ�С��λ����Ӧ��Ϊ������ݵ�λ������λ����һ��
% lb = min(data) �����½�
% ub = max(data) �����Ͻ�

%% ��������
if lb*ub <= 0 % �����Խ������ԭ��
    flag = 1; % ���ݻ�ͼ�淶��������ʾ��Χ�������ԭ��
else
    flag = 0;
end
range = ub-lb; % ���ݵļ���
maxData = max(abs([lb,ub])); % ����ֵ��������

%% �����Ϊ����
stringData = num2str(maxData); % ���������maxData��Ϊ�ַ���
if ~isempty(find(stringData=='.')) % �ж�maxData�Ƿ�ΪС����С����Ȼ����С����
    decimal_digit = length(stringData) - find(stringData=='.'); % С�����ֵ�λ��
    if stringData(1) == '0' % ��С����pure decimal��
        all_digit = length(stringData)-2; % ��Ϊ�������ȫ��λ��
    else % ��С����mixed decimal��
        all_digit = length(stringData)-1; % ��Ϊ�������ȫ��λ��
    end
    amp = decimal_digit;
    temp_range = 10^amp*range; % �����Ϊ����
    temp_digit = floor(log10(temp_range))+1; % �¼����λ��
    % ����������ʵ��N������Ϊmλ�������У�10^(m-1)<=N<10^m
    % �Ӷ���:m-1<=log10(N)<m �� log10(N)<m<=log10(N)+1��mΪ���� 
    % ��ˣ�mΪlog10(N)+1������ȡ��
    if temp_digit == 1
        amp = amp+1; % ���ڸ�λ������Ҫ��Ϊʮλ��
    end
else % maxDataΪ����
    digit = floor(log10(range))+1; % �����λ��
    all_digit = length(num2str(maxData)); % ��Ϊ�������ȫ��λ��
    if digit == 1
        amp = 1; % ���ڸ�λ������Ҫ��Ϊʮλ��
    else
        amp = 0;
    end
end
temp_range = round(10^amp*range); % ���ھ������⣬����ʹ��round
% ÿ�ݵ�λ��ֻ��Ϊm��m-1��
% ��-20~20���Ա���Ϊ-20:5:20��Ҳ���Ա���Ϊ-20:10:20
% ��ˣ���Ҫ������ʱ����10^amp���ٽ��������ж�

%% ÿ�ݵĴ�С��Ҫ��divisor_delta����
if all_digit >= 3
    amp1 = all_digit-2;
    divisor_delta = 5*10^amp1; 
else
    divisor_delta = 1;
end

%% ȷ��������ʾ��Χ
axis_display_data = [];
for divisor = [4,5,6] % �������ݷ�Ϊ�ķ���
    if rem(temp_range,divisor) == 0 % �����ܹ��ȷ�Ϊ���ɷ�
        temp_axis_display_data = linspace(lb,ub,divisor+1); % �ⶨ��������ʾ��Χ
        delta = temp_range/divisor; % ÿ�ݵĴ�С
        if rem(delta,divisor_delta) == 0 % ÿ�ݵĴ�С�ܹ���divisor_delta����
            delta = delta/10^amp; % delta����ʵֵ
            if flag == 1 % ������ʾ��Χ�б�Ҫ����ԭ��
                if ~isempty(find(temp_axis_display_data==0)) % ��ǰ�ⶨ��Χ����ԭ��
                    axis_display_data = temp_axis_display_data;
                    break
                end
            else % û�б�Ҫ����ԭ�� 
                axis_display_data = temp_axis_display_data;
                break
            end
        end
    end
end

if isempty(axis_display_data) % �ж��Ƿ���ȻΪ�ռ� 
    divisor = 4; % ����ȻΪ�ռ����ⶨ��Ϊ4��
    delta = fix(temp_range/divisor);
    temp = round(delta/divisor_delta)*divisor_delta; % ÿ�ݵĴ�С��Ҫ��divisor_delta����
    if temp ~= 0
        delta = round(delta/divisor_delta)*divisor_delta;
    else
        divisor_delta = divisor_delta/5;
        delta = round(delta/divisor_delta)*divisor_delta;
    end
    delta = delta/10^amp; % delta����ʵֵ
    if flag == 1 % ������ʾ��Χ�б�Ҫ����ԭ��
        num_ub = ceil(ub/delta); % ԭ�����ϵķ���
        num_lb = ceil(-lb/delta); % ԭ�����µķ���
        axis_display_data = -num_lb*delta:delta:num_ub*delta;
    else % û�б�Ҫ����ԭ��
        num = ceil((ub-lb)/delta); % ����
        lb = floor(10^amp*lb/divisor_delta)*divisor_delta/10^amp;
        if lb+num*delta >= ub
            axis_display_data = linspace(lb, lb+num*delta, num+1);
        else
            axis_display_data = linspace(lb, lb+(num+1)*delta, num+2);
        end
    end
end

%% ȷ��������ʾ����
stringData = num2str(delta); % ��ÿ�ݵĴ�Сdelta��Ϊ�ַ���
if ~isempty(find(stringData=='.')) % �ж�delta�Ƿ�ΪС����С����Ȼ����С����
    decimal_digit = length(stringData) - find(stringData=='.'); % С�����ֵ�λ��
    int_digit = find(stringData=='.')-1; % �������ֵ�λ��
else
    decimal_digit = 0;
    int_digit = length(stringData);
end