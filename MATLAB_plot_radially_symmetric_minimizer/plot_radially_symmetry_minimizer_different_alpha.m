clear all
close all
clc

epsilon = 0.0001; % a perturbation parameter to avoid NaN 

% specify parameters (require k*R > pi)
k = 1; 
a_array = linspace(0.5,4.5,6); 
R = 4; 

figure 
hold on 
for j=length(a_array):-1:1 
    a = a_array(j); 
    temp = @(R0) (-bessely(3/2,k*R0)*besselj(1/2,k*R) + besselj(3/2,k*R0)*bessely(1/2,k*R))/besselk(1/2,a*R) - (k/a)*(-bessely(3/2,k*R0)*besselj(3/2,k*R) + besselj(3/2,k*R0)*bessely(3/2,k*R))/besselk(3/2,a*R); 
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
    f = @(t) (t<R0) + (t>R0).*(t<R).*( ((k*R0)*sin(k*R0)+cos(k*R0))*sin(k*t)./(k*t) - (sin(k*R0) - (k*R0)*cos(k*R0))*cos(k*t)./(k*t) ) + (t>R).*( ((k*R0)*sin(k*R0)+cos(k*R0))*sin(k*R)./(k*R) - (sin(k*R0) - (k*R0)*cos(k*R0))*cos(k*R)./(k*R)   )*( a*R*exp(a*R) )*(exp(-a*t)/(a*t));
    fplot(f,[0,10])
end
grid on 
legend(['\alpha = ', num2str(a_array(6))],['\alpha = ', num2str(a_array(5))],['\alpha = ', num2str(a_array(4))],['\alpha = ', num2str(a_array(3))],['\alpha = ', num2str(a_array(2))],['\alpha = ', num2str(a_array(1))],'FontSize',16)
xlabel('|x|', 'FontSize',16)
ylabel('u_{*}^{rad}(x)','FontSize',16)
title(['d = 3, R = ', num2str(R), ', \kappa = ', num2str(k)],'FontSize',16)