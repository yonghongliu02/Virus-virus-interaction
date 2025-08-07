function ydot = f4(t,y,theta,xdata)

% INITIAL CONITIONS
N=2186*10000;

%virus1 cov; virus2 flu.

IS = theta(1);%只感染新冠
PS = 0; 
RS = theta(2);%新冠痊愈


SI = theta(3);%只感染流感
II = 0;
PI = 0;
RI = 0;

SP = 0;
IP = 0;
PP = 0;
RP = 0;

SR = theta(4);%流感痊愈
IR = 0;
PR = 0;
RR = theta(5);%都痊愈

SS = N - SI - IS - SR - RS - RR; 
I1=IS+II+IR+IP;
I2=SI+II+PI+RI;

%% Parameter

%beta
B1=theta(6);
B2=theta(7);

c1=theta(8);
c2=theta(9);

d1=theta(10);
d2=theta(11);


%sigma: interaction strength
sigma1=theta(12);%新冠对流感;
sigma2=theta(13);%流感对新冠

%p: temporary immunity period/
p1=1/theta(14);%新冠对流感的交叉保护时长
p2=1/theta(15);%流感对新冠的交叉保护时长


%gamma: recovery period
gamma1=1/5;%[PLoS Pathog. 2023 Mar 8;19(3):e1011167.]
gamma2=1/5;%【ref】

%w:Rate of waning immunity
w1=1/200;%新冠的保护时长
w2=1/365/2;%流感的保护时长

%% model

for i =2:180

%beta
beta1=B1*(1+c1*cos(2*pi*i/(365/2)-d1));
beta2=B2*(1+c2*cos(2*pi*i/365-d2));


% Susceptible to virus 1(cov), S (1st column)
SS(i) = SS(i-1) + w2*SR(i-1) + w1*RS(i-1) - beta1*SS(i-1)*I1(i-1)/N - beta2*SS(i-1)*I2(i-1)/N;
IS(i) = IS(i-1) + beta1*SS(i-1)*I1(i-1)/N + w2*IR(i-1) - beta2*(1-sigma1)*IS(i-1)*I2(i-1)/N - gamma1*IS(i-1);
PS(i) = PS(i-1) + gamma1*IS(i-1) + w2*PR(i-1) - p1*PS(i-1) -beta2*(1-sigma1)*PS(i-1)*I2(i-1)/N;
RS(i) = RS(i-1) + p1*PS(i-1) + w2*RR(i-1) - w1*RS(i-1) - beta2*RS(i-1)*I2(i-1)/N;



% Infected with virus 1, I(2r colum)
SI(i) = SI(i-1) + beta2*SS(i-1)*I2(i-1)/N + w1*RI(i-1) - beta1*(1-sigma2)*SI(i-1)*I1(i-1)/N - gamma2*SI(i-1);
II(i) = II(i-1) + beta2*(1-sigma1)*IS(i-1)*I2(i-1)/N + beta1*(1-sigma2)*SI(i-1)*I1(i-1)/N - gamma1*II(i-1) - gamma2*II(i-1);
PI(i) = PI(i-1) + gamma1*II(i-1) + beta2*(1-sigma1)*PS(i-1)*I2(i-1)/N  - p1*PI(i-1) - gamma2*PI(i-1);
RI(i) = RI(i-1) + p1*PI(i-1) +  beta2*RS(i-1)*I2(i-1)/N  - gamma2*RI(i-1)  - w1*RI(i-1);

%Post-infecte with virus 1 (3th column)
SP(i) = SP(i-1) + gamma2*SI(i-1) + w1*RP(i-1) - p2*SP(i-1) - beta1*(1-sigma2)*SP(i-1)*I1(i-1)/N;
IP(i) = IP(i-1) + beta1*(1-sigma2)*SP(i-1)*I1(i-1)/N + gamma2*II(i-1)  - gamma1*IP(i-1) - p2*IP(i-1);
PP(i) = PP(i-1) + gamma1*IP(i-1)  + gamma2*PI(i-1) - p1*PP(i-1)  - p2*PP(i-1);
RP(i) = RP(i-1) + p1*PP(i-1) +  gamma2*RI(i-1)  - w1*RP(i-1) - p2*RP(i-1);

% Recovere with virus 1 (4th column)
SR(i) = SR(i-1) + p2*SP(i-1)  + w1*RR(i-1) - beta1*SR(i-1)*I1(i-1)/N - w2*SR(i-1);
IR(i) = IR(i-1) + beta1*SR(i-1)*I1(i-1)/N + p2*IP(i-1)  - gamma1*IR(i-1) - w2*IR(i-1);
PR(i) = PR(i-1) + gamma1*IR(i-1) + p2*PP(i-1) - p1*PR(i-1) - w2*PR(i-1);
RR(i) = RR(i-1) + p1*PR(i-1) + p2*RP(i-1)  - w1*RR(i-1) - w2*RR(i-1);

%Infected with virus 1 or 2
%新冠
I1(i) = beta1*SS(i-1)*I1(i-1)/N + w2*IR(i-1)+...%IS
        beta2*(1-sigma1)*IS(i-1)*I2(i-1)/N + beta1*(1-sigma2)*SI(i-1)*I1(i-1)/N + ...%II
        beta1*(1-sigma2)*SP(i-1)*I1(i-1)/N + gamma2*II(i-1) + ...%IP
        beta1*SR(i-1)*I1(i-1)/N + p2*IP(i-1) ;%IR

%流感
I2(i) = beta2*SS(i-1)*I2(i-1)/N + w1*RI(i-1)+ ...%SI
        beta2*(1-sigma1)*IS(i-1)*I2(i-1)/N + beta1*(1-sigma2)*SI(i-1)*I1(i-1)/N+ ...%II
        gamma1*II(i-1) + beta2*(1-sigma1)*PS(i-1)*I2(i-1)/N + ...%PI
        p1*PI(i-1) +  beta2*RS(i-1)*I2(i-1)/N  ;%RI



    
end
 
%ydot=[I1(:)];
ydot = [I2(:), I1(:)];  % I2:流感, I1:新冠

global last_day_values 

last_day = [SS(end), IS(end), PS(end), RS(end), ...
            SI(end), II(end), PI(end), RI(end), ...
            SP(end), IP(end), PP(end), RP(end), ...
            SR(end), IR(end), PR(end), RR(end)];

last_day_values = last_day;  % Store the values

end

