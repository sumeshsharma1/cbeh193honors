load('xno','xno') %Expected value gpc
load('finalvar','finalvar') %Variance gpc
treal = [0:1:1000]; %Different t for gpc 
tautemp = randn(1,1000).*12 + 120;
x0= [0 0 0 0 0];
options = odeset('NonNegative', 1);
omega = 1;
for i = 1:1:1000
    f1 = @(t,x) [((0.1*omega*(x(5).^2))./(x(5).^2 + 10*omega^2)) - 0.1*x(1) + 0.01*omega;(4/tautemp(i))*(x(1) - x(2));(4/tautemp(i))*(x(2) - x(3));(4/tautemp(i))*(x(3) - x(4));(4/tautemp(i))*(x(4) - x(5))]; 
    f2 = @(t,x) [((20*omega*(x(5).^2))./(x(5).^2 + 10*omega^2)) - 0.1*x(1) + 0.01*omega;(4/tautemp(i))*(x(1) - x(2));(4/tautemp(i))*(x(2) - x(3));(4/tautemp(i))*(x(3) - x(4));(4/tautemp(i))*(x(4) - x(5))]; 
    [t,xmc] = ode45(f1,[0:1:200],x0,options);
    x1 = [xmc(201,1) xmc(201,2) xmc(201,3) xmc(201,4) xmc(201,5)];
    [t1, xmc1] = ode45(f2,[201:1:1000],x1,options);
    xnew(i,:) = [xmc((1:end),1);xmc1((1:end),1)]';
end

tmc = [t;t1];
%tmc = [tmc(1:10:end)]';

var = var(xnew); %Variance monte carlo
realvar = [finalvar(1:10:end)];
save('var','var')
figure
plot(treal,xno(:,1))
hold on
plot(treal,xno(:,6))
plot(treal,xno(:,11))
xlabel('Time')
ylabel('Expected Value')
plot(treal,finalvar)
plot(tmc,var)
legend('Sample = 0','Sample = sqrt(3)','Sample = -sqrt(3)','gPC Variance', 'Monte Carlo Variance')
title('Hermite Expected Value for Each Sample & gPC + Monte Carlo Variance')

    