n_simulations = 100;
a = .225; m = 33; b = 0.23; c = 0.23; p0= 10; %reput initial conditions here
tstep = .1;
t = 0:.1:1000;                     
steps = length(t);
x = zeros(10,steps,n_simulations);
xp = zeros(10,steps,n_simulations);
for i = 1:n_simulations
    x(1:10,1,i) = [10 1 0 0 0 0 0 0 0 0]';
    xp(1:10,1,i) = x(:,1,i);
    R = [];
    R = randn(1,2);
    %Make new parameter values with a distribution
    random_par = [.28 1.2].*R(1,2) + [2.8 12]; %Tp then Tm
    tau1 = random_par(1,1); %Tp
    tau2 = random_par(1,2); %Tm
    %Make new for loop for k = 1:timesteps
    for k = 2:steps
        %Make noise with zero mean and std sigma, multiplied by randn
        xp(1,k,i) = xp(1,k-1,i) + tstep*(a*x(10,k-1,i) - b*x(1,k-1,i)) + abs((sqrt(abs(a*x(10,k-1,i)))) - sqrt(abs(b*x(1,k-1,i))))*randn(1); %p(t)
        xp(2,k,i) = xp(2,k-1,i) + tstep*((m/(1 + ((x(6,k-1,i))^2)/(p0^2))) - c*x(2,k-1,i)) + abs((sqrt(abs((m/(1 + ((x(6,k-1,i))^2)/(p0^2))))))-sqrt(abs(c*x(2,k-1,i))))*randn(1); %m(t)
        xp(3,k,i) = xp(3,k-1,i) + tstep*((4/tau2)*(x(1,k-1,i) - x(3,k-1,i))); %p(t) time delay
        xp(4,k,i) = xp(4,k-1,i) + tstep*((4/tau2)*(x(3,k-1,i) - x(4,k-1,i)));
        xp(5,k,i) = xp(5,k-1,i) + tstep*((4/tau2)*(x(4,k-1,i) - x(5,k-1,i)));
        xp(6,k,i) = xp(6,k-1,i) + tstep*((4/tau2)*(x(5,k-1,i) - x(6,k-1,i)));
        xp(7,k,i) = xp(7,k-1,i) + tstep*((4/tau1)*(x(2,k-1,i) - x(7,k-1,i))); %m(t) time delay
        xp(8,k,i) = xp(8,k-1,i) + tstep*((4/tau1)*(x(7,k-1,i) - x(8,k-1,i)));
        xp(9,k,i) = xp(9,k-1,i) + tstep*((4/tau1)*(x(8,k-1,i) - x(9,k-1,i)));
        xp(10,k,i) = xp(10,k-1,i) + tstep*((4/tau1)*(x(9,k-1,i) - x(10,k-1,i)));
        %x1(k) = x1(k-1) + tstep(f1 + sqrt(abs(a*x(10,k-1) - (b*x(1,k-1))))*randn(0,1))
        x(:,k,i) = xp(:,k,i);
    end
end

x1sim = squeeze(xp(1,:,:));
x2sim = squeeze(xp(2,:,:));

var1 = var(x1sim');
var2 = var(x2sim');


figure
hold on; plot(t,var1); plot(t,var2)
legend('P(t)','M(t)')
% Make error plots, take variance vector for MC validation and EKF and subtract