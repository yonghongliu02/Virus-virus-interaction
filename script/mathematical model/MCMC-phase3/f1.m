clear model data params options
clear
filename = '/Users/lyh/Documents/CDC-QQB/合作单位资料/南安普顿大学/全球流感和新冠流行规律/Matlab/1-乘法模型估算/estimated_cases0731.csv';       %文件名
newdf=readtable(filename); 

%仓室值
comp=readtable('/Users/lyh/Documents/CDC-QQB/合作单位资料/南安普顿大学/全球流感和新冠流行规律/Matlab/phase3_86_360.xlsx'); 
%rowNames = comp.Properties.VariableNames';        % 原变量名作为行名
%comp=table2array(comp)';
%comp = array2table(comp, 'VariableNames', "comp", 'RowNames', rowNames);




t1=86+360+1;
t2=86+540;

newdf=newdf(t1:t2,:);


xdata.comp=comp;

data.xdata=xdata;
data.ydata = table2array(newdf(:,[5,2]));  % 4流感, 2新冠




%%
% The model sum of squares in file <algaess.html |algaess.m|> is
% given in the model structure.
model.ssfun = @f2;

%%
% All parameters are constrained to be positive. The initial
% concentrations are also unknown and are treated as extra parameters.

params1 = {
        %initial, min, max, pri_mu, pri_sig
       
     
        %参数beta：B1*(1+c1*cos(2*pi/365/2*(i-d1)));
        {'B1',   1.47,0,5,1.47,0.5}% B1=theta(6);
        {'B2',  0.5,0.25,2,0.5,0.5}% B2=theta(7);
        {'c1',   0.2,0,1,0.2,0.1}% c1=theta(8);
        {'c2',   0.3,0,1,0.3,0.1}% c2=theta(9);
        {'d1',   0,-60,60,0,1}% d1=theta(10);
        {'d2',   -4,-60,60,-4,1}% d2=theta(11);



        %sigma: interaction strength
        {'sigma1',   0,-0.05,0.05,0,0.2}% 新冠对流感;
        {'sigma2',   0.5,0,1,0.5,0.2}% 流感对新冠
        %p: temporary immunity period/
       {'p1',   1,1,10,1,2}% 新冠;
        {'p2',   30,1,60,30,2}% 流感;




        };
%%
% We assume having at least some prior information on the
% repeatability of the observation and assign rather non informational
% prior for the residual variances of the observed states. The default
% prior distribution is sigma2 ~ invchisq(S20,N0), the inverse chi
% squared distribution (see for example Gelman et al.). The 3
% components (_A_, _Z_, _P_) all have separate variances.

model.S20 = [4];
model.N0  = [1];

%%
% First generate an initial chain.
options.nsimu = 100000;
options.stats = 1;
[results, chain, s2chain]= mcmcrun(model,data,params1,options);


% % %%generate  chain after burn-in
options.nsimu = 100000;
options.stats = 1;
[results2, chain2, s2chain2] = mcmcrun(model,data,params1,options,results);

%%
% Chain plots should reveal that the chain has converged and we can
% % use the results for estimation and predictive inference.
figure
mcmcplot(chain2,[],results2); %,'pairs'
figure
mcmcplot(chain2,[],results2,'denspanel',2);


%%
% Function |chainstats| calculates mean ans std from the chain and
% estimates the Monte Carlo error of the estimates. Number |tau| is
% the integrated autocorrelation time and |geweke| is a simple test
% for a null hypothesis that the chain has converged.

results2.sstype = 1; % needed for mcmcpred and sqrt transformation

stats=chainstats(chain2,results2)


%%
% In order to use the |mcmcpred| function we need
% function |modelfun| with input arguments given as
% |modelfun(xdata,theta)|. We construct this as an anonymous function.

modelfun = @(d,th) f3(d(:,1),th,th(end),d);


% We sample 1000 parameter realizations from |chain| and |s2chain|
% and calculate the predictive plots.
nsample = 1000;
results2.sstype = 1;
out = mcmcpred(results2,chain2,s2chain2,data.xdata,modelfun,nsample);%data.ydata-->data



%% plot time series
% 流感拟合结果
tt = datetime(2022,11,7) + days((t1-1):(t2-1));
length=180;

figure
subplot(2,1,1)
fillyy(tt(1:length),out.obslims{1,1}{1,1}(3,:),out.obslims{1,1}{1,1}(1,:),[0.8 0.8 0.8]);
hold on 
plot(data.ydata(:,1),'.k');
plot(out.obslims{1,1}{1,1}(2,:));
hold off
title('流感拟合结果');
ylabel('病例数')

% 新冠拟合结果
subplot(2,1,2);
fillyy(tt(1:length), out.obslims{1,1}{1,2}(3,:), out.obslims{1,1}{1,2}(1,:), [0.8 0.8 0.8]);
hold on;
plot(data.ydata(:,2), '.k');  % 
plot(out.obslims{1,1}{1,2}(2,:)); 
title('新冠拟合结果');
ylabel('病例数');

%% output fitting parameter
CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));

Names = results2.names;
Mean = transpose(mean(chain2,1));
tmp=transpose(CIFcn(chain2,95));
UpCI = tmp(:,2);
DownCI = tmp(:,1);
geweke = stats(:,5);
T = table(Names,DownCI,Mean,UpCI, geweke, 'VariableNames', {'Parameter' ,'DownCI', 'Mean','UpCI','geweke'} );

%writetable(T, strcat(".\result\", cityname,'_fitting_parameter.csv'));

results2.dic


%% 输出16仓室
global last_day_values 

  % 定义各个仓室的名称
    compartment_names = {
        'SS', 'IS', 'PS', 'RS', ...
        'SI', 'II', 'PI', 'RI', ...
        'SP', 'IP', 'PP', 'RP', ...
        'SR', 'IR', 'PR', 'RR'
    };

    % 将 last_day 存储为表格
    last_day_table = array2table(last_day_values, 'VariableNames', compartment_names);

    writetable(last_day_table, '/Users/lyh/Documents/CDC-QQB/合作单位资料/南安普顿大学/全球流感和新冠流行规律/Matlab/phase4_86_540.xlsx');


%% 准备输出数据
% 创建时间序列
tt = datetime(2022,11,7) + days((t1-1):(t2-1));

% 流感拟合结果数据
flu_data = table(tt(1:length)', data.ydata(:,1), ...
    out.obslims{1,1}{1,1}(2,:)', ... % 拟合值
    out.obslims{1,1}{1,1}(1,:)', ... % 下限
    out.obslims{1,1}{1,1}(3,:)', ... % 上限
    'VariableNames', {'Date', 'Observed', 'Fitted', 'LowerCI', 'UpperCI'});

% 新冠拟合结果数据
covid_data = table(tt(1:length)', data.ydata(:,2), ...
    out.obslims{1,1}{1,2}(2,:)', ... % 拟合值
    out.obslims{1,1}{1,2}(1,:)', ... % 下限
    out.obslims{1,1}{1,2}(3,:)', ... % 上限
    'VariableNames', {'Date', 'Observed', 'Fitted', 'LowerCI', 'UpperCI'});



% 写入Excel文件
output_filename = '/Users/lyh/Documents/CDC-QQB/合作单位资料/南安普顿大学/全球流感和新冠流行规律/Matlab/fitting_phase4_540.xlsx';

% 写入流感拟合结果
writetable(flu_data, output_filename, 'Sheet', '流感拟合结果');

% 写入新冠拟合结果
writetable(covid_data, output_filename, 'Sheet', '新冠拟合结果');

% 写入参数估计结果
writetable(T, output_filename, 'Sheet', '参数估计结果');
