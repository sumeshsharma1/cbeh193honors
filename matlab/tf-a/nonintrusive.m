load('phiprod')
omega = 1; %can be varied
lambda = [1 0 -1; 1 sqrt(3) 2; 1 -sqrt(3) 2]; %roots of third degree Hermite+nominal solution
pi = (((transpose(lambda))*lambda)^-1)*transpose(lambda); %creation of pi matrix
tau1 = 1.176e-4*0 + (1/120);
tau2 = 1.176e-4*sqrt(3) + (1/120);
tau3 = 1.176e-4*-sqrt(3) + (1/120);
f1 = @(t,x) [((0.1*omega*(x(5).^2))./(x(5).^2 + 10*omega^2)) - 0.1*x(1) + 0.01*omega;(4*tau1)*(x(1) - x(2));(4*tau1)*(x(2) - x(3));(4*tau1)*(x(3) - x(4));(4*tau1)*(x(4) - x(5))]; 
f2 = @(t,x) [((0.1*omega*(x(5).^2))./(x(5).^2 + 10*omega^2)) - 0.1*x(1) + 0.01*omega; (4*tau2)*(x(1) - x(2)); (4*tau2)*(x(2) - x(3)); (4*tau2)*(x(3) - x(4)); (4*tau2)*(x(4) - x(5))];
f3 = @(t,x) [((0.1*omega*(x(5).^2))./(x(5).^2 + 10*omega^2)) - 0.1*x(1) + 0.01*omega; (4*tau3)*(x(1) - x(2)); (4*tau3)*(x(2) - x(3)); (4*tau3)*(x(3) - x(4)); (4*tau3)*(x(4) - x(5))];

% Construct functions for after t = 200
% Transform these
f4 = @(t,x) [((20*omega*(x(5).^2))./(x(5).^2 + 10*omega^2)) - 0.1*x(1) + 0.01*omega;(4*tau1)*(x(1) - x(2));(4*tau1)*(x(2) - x(3));(4*tau1)*(x(3) - x(4));(4*tau1)*(x(4) - x(5))]; 
f5 = @(t,x) [((20*omega*(x(5).^2))./(x(5).^2 + 10*omega^2)) - 0.1*x(1) + 0.01*omega; (4*tau2)*(x(1) - x(2)); (4*tau2)*(x(2) - x(3)); (4*tau2)*(x(3) - x(4)); (4*tau2)*(x(4) - x(5))];
f6 = @(t,x) [((20*omega*(x(5).^2))./(x(5).^2 + 10*omega^2)) - 0.1*x(1) + 0.01*omega; (4*tau3)*(x(1) - x(2)); (4*tau3)*(x(2) - x(3)); (4*tau3)*(x(3) - x(4)); (4*tau3)*(x(4) - x(5))];
%Solve before t = 200
x0 = [0 0 0 0 0]; %x0(1) is state, rest are time delay terms
options = odeset('NonNegative', 1);
[t, xno1] = ode45(f1,[0:1:200],x0,options);
[t, xno2] = ode45(f2,[0:1:200],x0,options);
[t, xno3] = ode45(f3,[0:1:200],x0,options);

%Construct initial condition vectors for after t = 200

x1 = [xno1(201,1) xno1(201,2) xno1(201,3) xno1(201,4) xno1(201,5)];
x2 = [xno2(201,1) xno2(201,2) xno2(201,3) xno2(201,4) xno2(201,5)];
x3 = [xno3(201,1) xno3(201,2) xno3(201,3) xno3(201,4) xno3(201,5)];
%Solve after t = 200
[t1, xno4] = ode45(f4,[201:1:1000],x1,options);
[t1, xno5] = ode45(f5,[201:1:1000],x2,options);
[t1, xno6] = ode45(f6,[201:1:1000],x3,options);
%Concatenate t < 200 and t > 200 matrices together
xno = [[xno1;xno4] [xno2;xno5] [xno3;xno6]];
%Concatenate time matrices
treal = [t;t1];
%Multiply pi by solution matrix
xeval = pi*[xno(:,1) xno(:,6) xno(:,11)]';
xnonew = [xno(:,1) xno(:,6) xno(:,11)]'; %for final var

mu = pi(1,:)*xnonew;
save('mu.mat','mu')

tmpreal = zeros(1,1001);
for i = 1:1:3
    for j = 1:1:3
        for k = 1:1:3
           tmp = pi(k,i)*pi(k,j)*phiprod(k)*xnonew(i,:).*xnonew(j,:);
           tmpreal = tmpreal+tmp;
        end
    end
end

finalvar = sum(tmpreal,1) - mu.^2;

figure
subplot(1,2,1),
plot(treal,xeval(1,:)) %plots non-intrusive solution
xlabel('Time')
ylabel('Expected Value')
title('Expected Value over Time: Hermite')
subplot(1,2,2),
plot(treal,finalvar); xlabel('Time'); ylabel('Variance'); title('Variance Over Time: Hermite')


save('xeval','xeval')
save('finalvar','finalvar')
save('xno','xno')
save('pi','pi')
save('xnonew','xnonew')



