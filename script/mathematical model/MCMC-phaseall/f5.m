%predict 

modelfun = @(d,th) f3(d(:,1),th,th(end),d);
nsample = 500;
results2.sstype = 1;
out1 = mcmcpred(results2,chain2,s2chain2,data.xdata,modelfun,nsample);
out2 = mcmcpred2(results2,chain2,s2chain2,data.xdata,modelfun,nsample);
%figure
%mcmcpredplot(out1,out.data,data);


%%
%出图
t1=87;
t2=86+540;

tt = datetime(2022,11,7) + days((t1-1):(t2-1));
length=540;

figure
subplot(2,1,1)
fillyy(tt(1:length),out1.obslims{1,1}{1,1}(3,:),out1.obslims{1,1}{1,1}(1,:),[0.8 0.8 0.8]);
hold on 
plot(data.ydata(:,1),'.k');
plot(out1.obslims{1,1}{1,1}(2,:));
hold off
title('流感拟合结果');
ylabel('病例数')

% 新冠拟合结果
subplot(2,1,2);
fillyy(tt(1:length), out1.obslims{1,1}{1,2}(3,:), out1.obslims{1,1}{1,2}(1,:), [0.8 0.8 0.8]);
hold on;
plot(data.ydata(:,2), '.k');  % 
plot(out1.obslims{1,1}{1,2}(2,:)); 
title('新冠拟合结果');
ylabel('病例数');


%% 准备输出数据
% 创建时间序列
tt = datetime(2022,11,7) + days((t1-1):(t2-1));

% 流感拟合结果数据
flu_data = table(tt(1:length)', data.ydata(:,1), ...
    out1.obslims{1,1}{1,1}(2,:)', ... % 拟合值
    out1.obslims{1,1}{1,1}(1,:)', ... % 下限
    out1.obslims{1,1}{1,1}(3,:)', ... % 上限
    'VariableNames', {'Date', 'Observed', 'Fitted', 'LowerCI', 'UpperCI'});

% 新冠拟合结果数据
covid_data = table(tt(1:length)', data.ydata(:,2), ...
    out1.obslims{1,1}{1,2}(2,:)', ... % 拟合值
    out1.obslims{1,1}{1,2}(1,:)', ... % 下限
    out1.obslims{1,1}{1,2}(3,:)', ... % 上限
    'VariableNames', {'Date', 'Observed', 'Fitted', 'LowerCI', 'UpperCI'});



% 写入Excel文件
output_filename = '/Users/lyh/Documents/CDC-QQB/合作单位资料/南安普顿大学/全球流感和新冠流行规律/Matlab/fitting_phase_all_simulation.xlsx';

% 写入流感拟合结果
writetable(flu_data, output_filename, 'Sheet', '流感拟合结果');

% 写入新冠拟合结果
writetable(covid_data, output_filename, 'Sheet', '新冠拟合结果');
