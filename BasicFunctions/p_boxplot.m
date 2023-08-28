function p_boxplot(data,lb,ub,names,subtitle,label_vertical)
columns = size(data,2);
figure;
boxplot(data,'labels',names);
number = length(names);
axis([0.5 number + 0.5 lb ub])
ax = gca;
ax.YGrid = 'on';
ylim1=get(gca,'ylim');
% set(gca,'YTick',93:0.5:100)
ylabel(label_vertical)
n_outliers = sum(data < lb);
for i = 1 : columns
    temp_outliers = n_outliers(i);
    if temp_outliers ~= 0
    str_temp = ['+', int2str(temp_outliers)];
    yvalue = lb + diff(ylim1)*.05;
    text(i + 0.1,yvalue,str_temp);
    end
end
title(subtitle)

x = 1:6;
colors = rand(6,3);
h = findobj(gca,'Tag','Box');
for j = 1:length(h)
    patch(get(h(j),'Xdata'), get(h(j),'Ydata'), colors(j,:), 'FaceAlpha',.5);
end