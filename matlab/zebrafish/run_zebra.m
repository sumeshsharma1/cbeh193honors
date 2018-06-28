t = [0:1000];
x0 = zeros(1,60);
x0(1) = 10; %p(t) = 10 at t = 0
x0(7) = 1; %m(t) = 1 at t = 0
options = odeset('NonNegative',1);
[t, xvalues] = ode15s('zebra_paper',t,x0);

plot(t,xvalues(:,1)) %plots p(t)
figure
plot(t,xvalues(:,7)) %plots m(t)