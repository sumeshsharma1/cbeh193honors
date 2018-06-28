k3 = .1066; %min^-1
k4 = .10658; %min^-1
tau = 6.4; %min
epo = 100;
x0 = [1 0 0 0 0 0 0 0]';
options = odeset('NonNegative',1);
n_simulations = 100;
k1 = randn(1,n_simulations)*.0021 + .021;
k2 = randn(1,n_simulations)*.246 + 2.46;
%k1 = .021*ones(1,n_simulations);
%k2 = 2.46*ones(1,n_simulations);
dt = .01;
t = 0:.01:10; % time vector
%g = .01;
steps = length(t);
z = [];
for i = 1:n_simulations
    f = @(t,x)[-k1(1,i)*x(1)*epo + 2*k4*x(8); -k2(1,i)*((x(2))^2) + k1(1,i)*x(1)*epo; -k3*x(3) + .5*k2(1,i)*((x(2))^2); k3*x(3) - k4*x(8); (4/tau)*(x(3)-x(5)); (4/tau)*(x(5)-x(6)); (4/tau)*(x(6) - x(7)); (4/tau)*(x(7)-x(8))];
    [t,xtemp] = ode15s(f,t,x0,options);
    %z(:,i,:) = xtemp;
end