clc
close all
clear all

% 1. Data description & Data plotting
load ("data.mat");
N = 100;
data = iddata(y, u);
nB = 40;

figure('Name', 'Input & Output');
subplot(2,1,1); 
plot(u); 
title('Input Signal u(t)');
subplot(2,1,2); 
plot(y); 
title('Output Signal y(t)');

% Construction of the phi matrix
col = [0; u(1:end-1)];
row = zeros(1,nB);
Phi = toeplitz(col, row);

% Cases
lambdas = [2, 500, 0.1, 2];
betas   = [0.75, 0.75, 0.75, 0.3];
cases   = length(lambdas);

% Types of kernel
kernel_types = {'TC', 'DI', 'SS'};
sigma0 = 0.58;
orders_arxReg = [0 nB 1];
t = 0:(nB-1);

for k_idx = 1:length(kernel_types)
    k_type = kernel_types{k_idx};
    
    % Create a figure for evry kind of kernel
    figure('Name', ['Kernel Type: ', k_type], 'Position', [100, 100, 1000, 800]);
    
    for c = 1:cases
        lambda = lambdas(c);
        beta = betas(c);
        
        K = zeros(nB, nB);
        for i = 1:nB
            for j = 1:nB
                if strcmp(k_type, 'TC')
                    K(i,j) = lambda * (beta^max(i,j));
                    
                elseif strcmp(k_type, 'DI')
                    if i == j
                        K(i,j) = lambda * (beta^i);
                    end
                    
                elseif strcmp(k_type, 'SS')
                    c_i = beta^i;
                    c_j = beta^j;
                    K(i,j) = lambda * ( (c_i * c_j * min(c_i,c_j))/2 - (min(c_i,c_j)^3)/6 );
                end
            end
        end
        
        K_inv = inv(K + 1e-9 * eye(nB)); 
        
        % Kernel-based estimation
        sigma2Kinv = (sigma0^2) * K_inv;
        Lambda_reg = norm(sigma2Kinv);
        R = sigma2Kinv / Lambda_reg;
        
        opt = arxOptions;
        opt.Regularization.Lambda = Lambda_reg;
        opt.Regularization.R = R; 
        
        m_arxReg = arx(data, orders_arxReg, opt);
        
        theta_est = m_arxReg.b(2:end)';
        
        % Calculates confidence intervals
        sigma2_hat = m_arxReg.NoiseVariance;
        P = K - K * Phi' * inv(Phi * K * Phi' + sigma2_hat * eye(N)) * Phi * K;
        std_theta = sqrt(abs(diag(P)));
        
        % Plotting
        subplot(2, 2, c);
        plot(t, theta0, 'g', 'LineWidth', 1.5); hold on;
        plot(t, theta_est, 'r', 'LineWidth', 1.5);
        
        % adding the 95% confidence
        plot(t, theta_est + 1.96*std_theta, 'r--');
        plot(t, theta_est - 1.96*std_theta, 'r--');
        
        title(sprintf('Case %d: \\lambda=%g, \\beta=%g', c, lambda, beta));
        legend('True \theta_0', 'Est m_{ARX}', '95% CI', 'Location', 'best');
        grid on;
        
    end
    sgtitle(sprintf('Performance of %s Kernel across different hyperparameters', k_type));
end