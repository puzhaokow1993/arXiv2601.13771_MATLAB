clear all
close all
clc

epsilon = 0.0001; % a perturbation parameter to avoid NaN 

% specify parameters (require k*R larger than first zero of J_0, which is approximately 2.4)
a = 5; 
k = 3; 

% lower bound -------------------------------------------------------------
R = 1; 
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
    LB = @(t) (t<R0) + (t>R0).*(t<R).*( -(pi/2)*k*R0.*bessely(1,k*R0).*besselj(0,k.*t) + (pi/2).*k.*R0.*besselj(1,k.*R0).*bessely(0,k.*t) ) + (t>R).*( -(pi./2).*k.*R0.*bessely(1,k*R0).*besselj(0,k.*R) + (pi./2).*k.*R0.*besselj(1,k*R0).*bessely(0,k*R)  ).*( besselk(0,a.*t)./besselk(0,a.*R) );

% upper bound -------------------------------------------------------------
R = 3; 
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
    UB = @(t) (t<R0) + (t>R0).*(t<R).*( -(pi./2).*k.*R0.*bessely(1,k.*R0).*besselj(0,k.*t) + (pi./2).*k.*R0.*besselj(1,k.*R0).*bessely(0,k.*t) ) + (t>R).*( -(pi/2).*k.*R0.*bessely(1,k.*R0)*besselj(0,k.*R) + (pi/2).*k.*R0.*besselj(1,k.*R0).*bessely(0,k.*R)  ).*( besselk(0,a.*t)./besselk(0,a.*R) );
   
% main loop ---------------------------------------------------------------

% Parameters
T= 5; % Truncation level 
h = 0.1; % FEM size 
tau = 0.1 * h^2;  % step size for stability
max_iter = 50; % maximum iteration 

% Initial guess 
u0 = @(x,y) UB(x.^2 + y.^2); 

grad_norms = zeros(max_iter,1);  % store gradient norms
u_prev = [];  % will hold solution from previous iteration

% Gradient-descent iterations 
for iter = 1:max_iter 

    % Create PDE model
    model = createpde();
    
    % Create geometry: unit circle
    g = @() circleg(T); % Helper function below
    geometryFromEdges(model, circleg(T)); 

    % Apply Dirichlet BC: u = 0 on boundary
    applyBoundaryCondition(model, 'dirichlet', ...
    'Edge', 1:model.Geometry.NumEdges, ...
    'u', 0);

    % specify coefficient 
    % q = @(location,state) (-(k.^2 + a.^2).*(location.x>-1).*(location.y>-1).*(location.x+location.y<2) + (a.^2)).*(state.u < 1) + 1/tau; % G = triangle 
    % q = @(location,state) (-(k.^2 + a.^2).*(location.x.^2 + (location.y/2).^2 < 1) + (a.^2)).*(state.u < 1) + 1/tau; % G = ellipse 
    % q = @(location,state) (-(k.^2 + a.^2).*(location.x>-1).*(location.x<1).*(location.y>-1).*(location.y<1) + (a.^2)).*(state.u < 1) + 1/tau; % G = square
    q = @(location,state) (-(k.^2 + a.^2).*(location.x.^2 + location.y.^2 < 4) + (k.^2 + a.^2).*(location.y.^2 + (location.x - 1.5).^2 < 0.25) + (a.^2)).*(state.u < 1) + 1/tau; % G = punctured ball 
    f = @(location,state) u0(location.x,location.y)./tau; 

    specifyCoefficients(model, 'm', 0, 'd', 0, 'c', 1, 'a', q, 'f', f); 

    % Generate mesh
    generateMesh(model, 'Hmax', h, 'GeometricOrder', 'linear'); % smaller Hmax -> finer mesh

    % Solve PDE
    result = solvepde(model);
    u = result.NodalSolution;  % vector of nodal values 
    
    % ----- compute gradient norm -----
    if iter > 1
        diff_u = u - u_prev;  % change between iterations
        grad_norms(iter) = norm(diff_u, 2) / tau;  % approximate gradient norm
    end
    u_prev = u;

    % Interpolate the solution for next iteration 
    nodes = model.Mesh.Nodes;     % 2 × N array: [x; y]
    x_nodes = nodes(1, :)';
    y_nodes = nodes(2, :)';
    u_nodes = u(:);               % Ensure it's a column vector

    F = scatteredInterpolant(x_nodes, y_nodes, u_nodes, 'linear', 'none');
    u0 = @(x, y) max(min(F(x, y),UB(x.^2 + y.^2)),LB(x.^2 + y.^2)); % overwrite u0 for next iteration 

end

% Define a grid over your domain
x_vals = linspace(-2, 2, 100);  
y_vals = linspace(-2, 2, 100);
[X, Y] = meshgrid(x_vals, y_vals);

% Evaluate the function handle
Z = u0(X, Y);  % If function returns NaN outside domain, mask it later

% Plot solution 
figure;
surf(X, Y, Z, 'EdgeColor', 'none');
xlabel('x'), ylabel('y'), zlabel('u(x, y)');
colormap(jet); colorbar; shading interp; 
view(2);  % Optional: top view

% Plot gradient norm over iterations
figure;
semilogy(1:max_iter, grad_norms, '-o');
xlabel('Iteration');
ylabel('Approximate gradient norm');
grid on;