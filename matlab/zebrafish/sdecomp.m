%% ZebraFish SDE:
%Run more samples, 1000
%One w/o parametric uncertainties, one with 
n_species = 2;
n_simulations = 1000;

R = zeros(n_simulations,2);

par = [2.8;12]; % Nominal value, 2x1, taus
par_lb = [2.52;10.8]'; %Tp first row, Tm second row
par_ub = [3.08;13.2]'; %tau upper bound

x0 = [10;1;0;0;0;0;0;0;0;0]; %p0, m0, time dealys?!??!? 2x1

%x0_ub = [110;110]'; 
%x0_lb = [90;90]';

%random_IC = zeros(n_simulations,2);
random_par = zeros(n_simulations,2);

a = .225; k = 33; b = 0.23; c = 0.23; p0= 10; %reput initial conditions here
dt = .01;
t = 0:.01:1000;                     % time vector
steps = length(t);

for i=1:n_simulations
    R = [];
    R = randn(1,2);
    %random_IC(i,1:2) = (x0_ub - x0_lb).*R(i,1:2) + x0_lb;
    %x0(2) = random_IC(i,1);
    %x0(5) = random_IC(i,2);
    %random_par(i,1:2) = (par_ub - par_lb).*R(i,1:2) + par_lb;
    random_par = [.28 1.2].*R(1,1:2) + [2.8 12]; %Tp then Tm
    %Tp = random_par(1,1);
    %Tm = random_par(1,2);
    Tp = 2.8;
    Tm = 12;
    f = @(t,x)[a*x(10) - 0.23*x(1); (33/(1 + (x(6).^2)/100)) - 0.23*x(2); (4/Tm)*(x(1)-x(3)); (4/Tm)*(x(3)-x(4)); (4/Tm)*(x(4)-x(5)); (4/Tm)*(x(5)-x(6)); (4/Tp)*(x(2)-x(7)); (4/Tp)*(x(7)-x(8)); (4/Tp)*(x(8)-x(9)); (4/Tp)*(x(9)-x(10))];
    %g = @(t,x)[(sqrt(abs(a*x(10))) - b*x(1));(sqrt(abs((k/(1 + ((x(6))^2)/(p0^2))) - c*x(10))));0;0;0;0;0;0;0;0];
%     g = @(t,x)[(sqrt(abs(a*x(10)))) - sqrt(abs(b*x(1)));(sqrt(abs((k/(1 + ((x(6))^2)/(p0^2))))))-sqrt(abs(c*x(2)));0;0;0;0;0;0;0;0];
     g = @(t,x)[abs((sqrt(abs(a*x(10)))) - sqrt(abs(b*x(1))));abs((sqrt(abs((k/(1 + ((x(6))^2)/(p0^2))))))-sqrt(abs(c*x(2))));0;0;0;0;0;0;0;0];

    opts = sdeset('RandSeed',i);  % Set random seed
    y = sde_euler(f,g,t,x0,opts); % Integrate
    trajectory(:,:,i) = y;
    %i
end


sample_1_sde = squeeze(trajectory(:,1,:)); %p(t)
sample_2_sde = squeeze(trajectory(:,2,:)); %m(t)


% %% Plot the results here: 
% 
% for i = 1:n_simulations
%     plot(t,trajectory(:,:,i));
%     hold on;
% end
% hold off;

%% Calculate the variance throughout the whole stochastic trajectory:

% E_value = sum(trajectory)./n_simulations;
% 
% for i=1:steps
%     for j = 1:n_species
%         E_value_sde(i,j) = sum(trajectory(i,j,:))./n_simulations;
%     end
% end
% % 
% 
% for i=1:steps
%     for j = 1:n_species
%         Var_sde(i,j)  = sum(((trajectory(i,j,:) - E_value_sde(i,j))).^2)./n_simulations;
%     end
% end