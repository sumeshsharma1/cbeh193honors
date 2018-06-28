load('C')
load('phieval')
load('phiprod')
x0 = zeros(1,48);
x0(1) = 1;
epo = 100;
options = odeset('NonNegative',1);
[t, xvalues] = ode15s('jakstat_c',[0:.01:5],x0,options,epo);

xfin1 = xvalues(:,1:6)*phieval;
xfin2 = xvalues(:,7:12)*phieval;
xfin3 = xvalues(:,13:18)*phieval;
xfin4 = xvalues(:,19:24)*phieval;

figure
plot(t,xvalues(:,1)); hold on; plot(t,xvalues(:,7)); plot(t,xvalues(:,13)); plot(t,xvalues(:,19))
xlabel('Time'); ylabel('Expected Value'); title('Expected Value via gPC')
legend('X1','X2','X3','X4')

varphi = [];

for i = 1:1:6 %computes variance column vectors 
    varphi(i) = (phiprod(i)).^2;
    x1var(:,i) = ((xvalues(:,i)).^2).*varphi(i);
    x2var(:,i) = ((xvalues(:,i+6)).^2).*varphi(i);
    x3var(:,i) = ((xvalues(:,i+12)).^2).*varphi(i);
    x4var(:,i) = ((xvalues(:,i+18)).^2).*varphi(i);
end

x1var = sum(x1var,2) - ((xvalues(:,1)).^2); %sums vectors to get overall variance, fix
x2var = sum(x2var,2) - ((xvalues(:,7)).^2);
x3var = sum(x3var,2) - ((xvalues(:,13)).^2);
x4var = sum(x4var,2) - ((xvalues(:,19)).^2);

figure 
plot(t,x1var); hold on; plot(t,x2var); plot(t,x3var); plot(t,x4var)
xlabel('Time'); ylabel('Variance'); title('Variance via gPC')
legend('X1','X2','X3','X4')

    


