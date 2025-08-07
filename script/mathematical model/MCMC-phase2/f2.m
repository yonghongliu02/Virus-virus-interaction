function ss = f2(theta, data)
    time = (1:180)';
    ydata = data.ydata; % 现在为两列 [流感, 新冠]
    xdata = data.xdata;
    y0 = theta(end);
    ymodel = f3(time, theta, y0, xdata); % ymodel应为两列
    
    % 计算两列残差的平方和:平方根变换
  % ss = sum( (1/mean(ydata(:,1))) * (sqrt(ymodel(:,1)) - sqrt(ydata(:,1))).^2 ) + ...
  %   sum( (1/mean(ydata(:,2))) * (sqrt(ymodel(:,2)) - sqrt(ydata(:,2))).^2 );

   % 计算两列残差的平方和:对数变换
    %ss = sum( (1/mean(ydata(:,1))) * (log(ymodel(:,1)+1) - log(ydata(:,1)+1)).^2 ) + ...
    % sum( (1/mean(ydata(:,2))) * (log(ymodel(:,2)+1) - log(ydata(:,2)+1)).^2 );

 % 计算两列残差的平方和:基本
     ss = sum((sqrt(ymodel(:,1)) - sqrt(ydata(:,1))).^2) + ...
         sum((sqrt(ymodel(:,2)) - sqrt(ydata(:,2))).^2);
end