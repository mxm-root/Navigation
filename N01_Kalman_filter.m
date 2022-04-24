%==========================================================================
% Tutorial Navigation
% Topic : Kalman filter
% Author: M.Trifonov 
% Email: trifonov.m@yahoo.com
% Date(dd-mm-yyyy): 14-09-2020
%==========================================================================
clc; clear; close all
% 1. Initial data
randn('state',sum(100*clock)); % randomizing
X0=[3000;500;10;-10]; % actual initial state (initial position)
dt=5; % step time
Fi=eye(4,4); Fi(1,3)=dt; Fi(2,4)=dt; % transition matrix
C=zeros(2,4); C(1,1)=1; C(2,2)=1; % measurement matrix
% Measurment noise characteristics 
D_eta=900; Sig_eta=sqrt(D_eta); % noise variance & standard deviation - Sigma
K_eta=eye(2,2); K_eta=D_eta*K_eta; 
D_V0=25; Kx0=zeros(4,4); Kx0(1,1)=D_eta; Kx0(2,2)=D_eta; % initial state covariance
Kx0(3,3)=D_V0;Kx0(4,4)=D_V0;
N=25; % number of steps
X=X0;
%  2. Recurrent simulation & processing
for i=1:N % cycle in steps (counter)
    t(i) = i*dt; % time
    X = Fi*X; % actual current state
    X_data(i,1) = X(1);
    X_data(i,2) = X(2);
    y  = C*X; % ideal measurement (without error)
    eta = Sig_eta*randn(2,1);% measurement error
    y_izm = y + eta; % actual measurement (with error)
    Y_data(i,1) = y_izm(1);
    Y_data(i,2) = y_izm(2);
    % Forecast phase
    if i==1
        Papr = Fi*Kx0*Fi';
        Xapr = [y_izm;0.1;0.1];
    else
        % Covariance matrix forecast 
        Papr = Fi*Paps*Fi'; % before next observation
        % Estimation forecast
        Xapr = Fi*Xaps;
    end
    % Re-calculation phase
    Xapr_data(i,1)=Xapr(1); 
    Xapr_data(i,2)=Xapr(2);
    Paps=Papr-Papr*C'*inv(K_eta+C*Papr*C')*C*Papr; % after first observation
    Sxaps(i)=sqrt(Paps(1,1));
    SVxaps(i)=sqrt(Paps(3,3));
    Xaps=Xapr+Paps*C'*inv(K_eta)*(y_izm-C*Xapr);
    Eps=X-Xaps; epsx(i)=Eps(1);epsVx(i)=Eps(3);
    Xaps_data(i,1)=Xaps(1);
    Xaps_data(i,2)=Xaps(2);
end
disp(Eps);
% Plotting
figure(1)
grid on;hold on;
plot(X0(1),X0(2),'oR'); % actual initial position
axis ij, axis([2800 3200 -100 600]),axis equal;
plot(X_data(:,1),X_data(:,2),'DG'); % actual current position
plot(Y_data(:,1),Y_data(:,2),'xB'); % measurement points
plot(Xapr_data(:,1),Xapr_data(:,2),'c*'); % a priori trajectory estimation
plot(Xaps_data(:,1),Xaps_data(:,2),'SK'); % a posteriori trajectory estimation
legend('actual initial state','actual trajectory','trajectory from measurements','a priori trajectory estimation','a posteriori trajectory estimation');
xlabel('\itx \rmposition (m)');ylabel('\itz \rmposition (m)');
figure(2)
grid on;hold on;
stem(epsx);
plot(3*Sxaps,'--r'), legend({'$_{\hat{X}}$','$\pm3\sigma_{\hat{X}}$'},'Interpreter','latex');
plot(-3*Sxaps,'--r');
xlabel('Number of measurments \itN'); ylabel('Error of \itx \rmposition \rm(m)');
figure(3)
grid on;hold on;
stem(epsVx);
plot(3*SVxaps,'--r');
legend({'$_{\hat{V_x}}$','$\pm3\sigma_{\hat{V_x}}$'},'Interpreter','latex');
plot(-3*SVxaps,'--r');
xlabel('Number of measurments \itN');ylabel('Error of speed \itV_x \rm(m/s)');