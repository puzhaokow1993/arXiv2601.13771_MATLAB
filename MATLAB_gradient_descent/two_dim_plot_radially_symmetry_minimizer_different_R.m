clear all
close all
clc

epsilon = 0.0001; % a perturbation parameter to avoid NaN 

% specify parameters (require k*R larger than first zero of J_0, which is approximately 2.4)
R_array = [1,3]; 
a = 5; 
k = 3; 

figure 
hold on 
for j=1:length(R_array) 
    R = R_array(j); 
    temp = @(R0) (-bessely(1,k*R0)*besselj(0,k*R) + besselj(1,k*R0)*bessely(0,k*R))/besselk(0,a*R) - (k/a)*(-bessely(1,k*R0)*besselj(1,k*R) + besselj(1,k*R0)*bessely(1,k*R))/besselk(1,a*R); 
    x = linspace(0.0001,R-0.0001,1000); 
    y = temp(x); 
    idx = find(diff(sign(y))); % Find intervals where sign changes
    zeros_list = zeros(size(idx));
    for i = 1:length(idx)
        left = x(idx(i));
        right = x(idx(i)+1);
        zeros_list(i) = fzero(temp, [left, right]);
    end 
    R0 = max(zeros_list); % Get the largest one
    f = @(t) (t<R0) + (t>R0).*(t<R).*( -(pi/2)*k*R0*bessely(1,k*R0)*besselj(0,k*t) + (pi/2)*k*R0*besselj(1,k*R0)*bessely(0,k*t) ) + (t>R).*( -(pi/2)*k*R0*bessely(1,k*R0)*besselj(0,k*R) + (pi/2)*k*R0*besselj(1,k*R0)*bessely(0,k*R)  )*( besselk(0,a*t)/besselk(0,a*R) );
    fplot(f,[0,5])
end
grid on 
legend(['R = ', num2str(R_array(1))],['R = ', num2str(R_array(2))],'FontSize',16)
xlabel('|x|', 'FontSize',16)
ylabel('u_{*}^{rad}(x)','FontSize',16)
title(['d = 2, \beta = ', num2str(k), ', \alpha = ', num2str(a)],'FontSize',16)