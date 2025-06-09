clear;
close all;
clc;
rng(124)

%% Param Declaration
dt      = 0.01;
T       = 0:dt:25;
StepNum = numel(T);

% Actual Values of the Parameters
ATrue     = [1.0 0.3; 0.4 0.8];
GammaTrue = 0.3;

% Measurement Noise
SigmaY    = 0.05;  

%% Noise-Free State Trajectories
% Mem Alloc
X                 = zeros(2,StepNum);
X(:, 1)           = [0.10; 0.04];       % Init Cond

for k = 2:StepNum
    xDot         = (1 - X(:, k-1)) .* (ATrue*X(:, k-1)) - GammaTrue*X(:, k-1);
    X(:, k)      = X(:, k-1) + dt*xDot;
end

YMeas            = X + SigmaY*randn(size(X));           % Additive Noise

%% EKF Init
% Mem Alloc
Xhat        = zeros(7, StepNum);      
Xhat(:, 1)  = [0.08; 0.02; 0.6; 0.6; 0.6; 0.6; 0.5];    % Init Cond

P           = diag([1e-3, 1e-3, 0.3, 0.3, 0.3, 0.3, 0.1]);

% Process Noise
QState      = diag([1e-4, 1e-4]);
QParams     = diag(5e-8*ones(5,1));
Q           = blkdiag(QState, QParams);

% Measurement Noises
R           = SigmaY^2 * eye(2);         % Measurement Covariance
H           = [eye(2), zeros(2, 5)];     % Measurement Jacobian  (2×7)

%% EKF Main Loop
for k = 2:StepNum
    XPrev            = Xhat(:, k-1);

    % ---- Prediction ----
    XPred            = XPrev  + dt * ODE(XPrev);
    F                = eye(7) + dt * Jacobian(XPrev);   % discrete-time F ≈ I+ΔtA
    PPred            = F*P*F' + dt * Q;

    % ---- Update ----
    MeasRes          = YMeas(:, k) - XPred(1:2);
    S                = H * PPred * H' + R;
    K                = PPred * H' / S;                  % Kalman Gain
    Xhat(:, k)       = XPred + K * MeasRes;
    P                = (eye(7) - K*H)*PPred;
end

%% Plot Results
% Noise Free Trajectories
figure
subplot(2, 1, 1)
plot(T, X(1, :), '--', T, Xhat(1, :), 'LineWidth', 1.2)
ylabel('x_1  (node 1)')
legend({'True', 'Estimate'}, 'Location', 'best')
subplot(2,1,2)
plot(T, X(2, :), '--', T, Xhat(2, :), 'LineWidth', 1.2)
ylabel('x_2  (node 2)')
xlabel('Time  [s]')
sgtitle('Infection Levels')

% Parameter Estimates with References
ParNames   = {'a_{11}','a_{12}','a_{21}','a_{22}','\gamma'};
TrueParams = [ATrue(1, 1), ATrue(1, 2), ATrue(2, 1), ATrue(2, 2), GammaTrue];
ParamIdx   = 3:7;

figure
for p = 1:5
    subplot(5, 1, p)
    plot(T, Xhat(ParamIdx(p),:), 'LineWidth', 1.2)
    hold on
    yline(TrueParams(p), '--');
    ylabel(ParNames{p});
    ylim([0, 1])
    
    if p == 5
        xlabel('Time [s]')
    end
end
sgtitle('Parameter Estimates vs. Truth');

% Error Histograms
paramErrors  = Xhat(ParamIdx, :) - TrueParams';
figure;
for p = 1:5
    subplot(3,2,p);
    histogram(paramErrors(p, :), 50, 'Normalization', 'pdf');
    title(['Error ', ParNames{p}]);
    xlabel('Estimate − Truth'); 
    ylabel('pdf');
end
sgtitle('Parameter Error Distributions');
