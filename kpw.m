%% time scale for natural kalman/kalman plus weight/ reduced kalman plus weight
%edit by shixiaohui@foxmail.com 20181220
%according to C. A. Greenhall: A review of reduced kalman filters for clock ensembles
clc;clear all;

q_all = [10 1e-5 1e-8 1 1e-3 1e-6 1 1e-3 1e-5 1e-3 1e-3 1e-4 1e-3 1e-3 1e-4]'*1e-10;
N_clock = length(q_all)/3;%number of clocks
N_time = 1e4;%simulation count
delta = 1;%time duration
phi_delta = [1,delta, delta^2/2; 0, 1, delta; 0, 0, 1];%传递函数
phi = kron(eye(N_clock), phi_delta);
Q = diag(q_all)*0.16;
wi_tmp = [delta, delta^3/3, delta^5/20]*reshape(q_all,3, N_clock);
wi = 1./sum(1./(wi_tmp))./(wi_tmp);

x_raw = zeros(3*N_clock, N_time);
x_p = x_raw;
x_kpw = zeros(1,N_time);

%观测矩阵H
H = zeros(N_clock-1,N_clock*3);
H(:,1) = 1;
for n = 1:N_clock-1
    H(n,3*n+1) = -1;
end
R = eye(N_clock-1)*0;%测量方差
P = zeros(size(Q));%协方差矩阵
flag_reduced_kpw = true;
for t = 2:N_time
    %演化
    x_raw(:,t) = phi*x_raw(:,t-1)+q_all.*(rand(3*N_clock,1)-0.5);
    
    %预测
    xp_tmp = phi*x_p(:,t-1);
    P_tmp = phi*P*(phi') + Q;
    Z_tmp = H*xp_tmp;%预测的测量结果
    Z = H*x_raw(:,t);%真实测量结果
    I_tp1 = Z-Z_tmp;%真实测量结果和预测测量结果的误差
    C = H*P_tmp*H'+R;
    K = P_tmp*(H')*pinv(C);
    
    %更新变量和状态
    x_p(:,t) = xp_tmp+K*I_tp1;
    P =  (eye(N_clock*3)-K*H)*P_tmp;
    if flag_reduced_kpw
        P(1:3:end,:) = 0;
        P(:,1:3:end) = 0;
    end
    x_kpw(t) = x_kpw(t-1)+sum(wi*(x_raw(1:3:end,t)-x_raw(1:3:end,t-1)-delta*x_p(2:3:end,t-1)-1/2*delta^2*x_p(3:3:end,t-1)));
end
x_nk = x_raw(1,:) - x_p(1,:);
% x_nk = x(4,:) - x_p(4,:);%the same with x(1,:) - x_p(1,:)

%% plot the time domain result
figure
time = 1:N_time;
plot(time, x_raw(1:3:end,:));
hold on
plot(time, x_nk,'p',time, x_kpw,'x');

%% plot the allan deviation
x_phase =  [x_raw(1:3:end,:);x_nk;x_kpw];
tau_in = [1,2,4, 8, 10, 20 40 80 100 200 400 800 1000 2e3 4e3 8e3 1e4 2e4 4e4 8e4]; 
for n = 1:N_clock+2
    data.phase = x_phase(n,:);
    data.rate = 1/delta;
    [ad(:,n),~,~,tau] = allan(data,tau_in,[],0);
end
figure;
loglog(tau,ad(:, 1:N_clock),'p');
hold on
loglog(tau,ad(:, N_clock+1:N_clock+2),'x-');
hold off
grid on
