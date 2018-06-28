k3 = .1066; %min^-1
k4 = .10658; %min^-1
tau = 6.4; %min
epo = 100;
x0 = [1 0 0 0 0 0 0 0]';
n_simulations = 1000;
k1 = randn(1,n_simulations)*.0021 + .021;
k2 = randn(1,n_simulations)*.246 + 2.46;
%k1 = .021*ones(1,n_simulations);
%k2 = 2.46*ones(1,n_simulations);
dt = .01;
t = 0:.01:5; % time vector
%g = .01;
steps = length(t);
for i = 1:n_simulations
    f = @(t,x)[-k1(1,i)*x(1)*epo + 2*k4*x(8); -k2(1,i)*((x(2))^2) + k1(1,i)*x(1)*epo; -k3*x(3) + .5*k2(1,i)*((x(2))^2); k3*x(3) - k4*x(8); (4/tau)*(x(3)-x(5)); (4/tau)*(x(5)-x(6)); (4/tau)*(x(6) - x(7)); (4/tau)*(x(7)-x(8))];
    g = @(t,x)[sqrt(.005)*sqrt(-k1(1,i)*x(1)*epo + 2*k4*x(8));sqrt(.005)*sqrt(-k2(1,i)*((x(2))^2) + k1(1,i)*x(1)*epo);sqrt(.005)*sqrt(-k3*x(3) + .5*k2(1,i)*((x(2))^2));sqrt(.005)*sqrt(k3*x(3) - k4*x(8));0;0;0;0];
    %g = @(t,x)[sqrt(.005)*sqrt(x(1));sqrt(.005)*sqrt(x(2));sqrt(.005)*sqrt(x(3));sqrt(.005)*sqrt(x(4));0;0;0;0];
    opts = sdeset('RandSeed',i);  % Set random seed
    y = sde_euler(f,g,t,x0,opts); % Integrate
    trajectory(:,:,i) = y;
end

sample_1_sde = squeeze(trajectory(:,1,:)); %x1 
sample_2_sde = squeeze(trajectory(:,2,:)); %x2
sample_3_sde = squeeze(trajectory(:,3,:)); %x3
sample_4_sde = squeeze(trajectory(:,4,:)); %x4

covmat = [];
for j = 1:1:steps
    c12 = cov(sample_1_sde(j,:),sample_2_sde(j,:));
    c13 = cov(sample_1_sde(j,:),sample_3_sde(j,:));
    c14 = cov(sample_1_sde(j,:),sample_4_sde(j,:));
    c23 = cov(sample_2_sde(j,:),sample_3_sde(j,:));
    c24 = cov(sample_2_sde(j,:),sample_4_sde(j,:));
    c34 = cov(sample_3_sde(j,:),sample_4_sde(j,:));
    covmat(:,:,j) = [var(sample_1_sde(j,:)),c12(1,2),c13(1,2),c14(1,2);
        c12(1,2),var(sample_2_sde(j,:)),c23(1,2),c24(1,2);
        c13(1,2),c23(1,2),var(sample_3_sde(j,:)),c34(1,2);
        c14(1,2),c24(1,2),c34(1,2),var(sample_4_sde(j,:))];
end


var1 = var(sample_1_sde');
var2 = var(sample_2_sde');
var3 = var(sample_3_sde');
var4 = var(sample_4_sde');
figure; hold on
plot(t,var1), plot(t,var2), plot(t,var3), plot(t,var4)
xlabel('Time'); ylabel('Variance')
title('1000 MC Simulations of SDE JAK/STAT Solver')
legend('X1','X2','X3','X4')
    