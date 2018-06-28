load('phiprod')
%assume Ns = 7
epo = 100;
xi1 = [-sqrt(3) -sqrt(3) -sqrt(3) 0 0 0 sqrt(3) sqrt(3) sqrt(3)];
xi2 = [-sqrt(3) 0 sqrt(3) -sqrt(3) 0 sqrt(3) -sqrt(3) 0 sqrt(3)];

lambda = [];
phi0 = ones(length(xi1),1);
for i = 1:length(xi1)
    phi1(i) = xi1(i);
    phi2(i) = xi2(i);
    phi3(i) = (xi1(i))^2 - 1;
    phi4(i) = xi1(i)*xi2(i);
    phi5(i) = (xi2(i))^2 - 1;
end
lambda = [phi0 phi1' phi2' phi3' phi4' phi5'];

pi = (((transpose(lambda))*lambda)^-1)*transpose(lambda);
save('pi','pi')
k1 = xi1*.0021 + .021;
k2 = xi2*.246 + 2.46;

f0 = @(t,x) [-k1(1)*x(1)*epo + 2*.10658*x(8); -k2(1)*(x(2))^2 + k1(1)*x(1)*epo; -.1066*x(3) + .5*k2(1)*((x(2))^2); .1066*x(3) - .10658*x(8); (4/6.4)*(x(3) - x(5)); (4/6.4)*(x(5) - x(6)); (4/6.4)*(x(6) - x(7)); (4/6.4)*(x(7) - x(8))];
f1 = @(t,x) [-k1(2)*x(1)*epo + 2*.10658*x(8); -k2(2)*(x(2))^2 + k1(2)*x(1)*epo; -.1066*x(3) + .5*k2(2)*((x(2))^2); .1066*x(3) - .10658*x(8); (4/6.4)*(x(3) - x(5)); (4/6.4)*(x(5) - x(6)); (4/6.4)*(x(6) - x(7)); (4/6.4)*(x(7) - x(8))];
f2 = @(t,x) [-k1(3)*x(1)*epo + 2*.10658*x(8); -k2(3)*(x(2))^2 + k1(3)*x(1)*epo; -.1066*x(3) + .5*k2(3)*((x(2))^2); .1066*x(3) - .10658*x(8); (4/6.4)*(x(3) - x(5)); (4/6.4)*(x(5) - x(6)); (4/6.4)*(x(6) - x(7)); (4/6.4)*(x(7) - x(8))];
f3 = @(t,x) [-k1(4)*x(1)*epo + 2*.10658*x(8); -k2(4)*(x(2))^2 + k1(4)*x(1)*epo; -.1066*x(3) + .5*k2(4)*((x(2))^2); .1066*x(3) - .10658*x(8); (4/6.4)*(x(3) - x(5)); (4/6.4)*(x(5) - x(6)); (4/6.4)*(x(6) - x(7)); (4/6.4)*(x(7) - x(8))];
f4 = @(t,x) [-k1(5)*x(1)*epo + 2*.10658*x(8); -k2(5)*(x(2))^2 + k1(5)*x(1)*epo; -.1066*x(3) + .5*k2(5)*((x(2))^2); .1066*x(3) - .10658*x(8); (4/6.4)*(x(3) - x(5)); (4/6.4)*(x(5) - x(6)); (4/6.4)*(x(6) - x(7)); (4/6.4)*(x(7) - x(8))];
f5 = @(t,x) [-k1(6)*x(1)*epo + 2*.10658*x(8); -k2(6)*(x(2))^2 + k1(6)*x(1)*epo; -.1066*x(3) + .5*k2(6)*((x(2))^2); .1066*x(3) - .10658*x(8); (4/6.4)*(x(3) - x(5)); (4/6.4)*(x(5) - x(6)); (4/6.4)*(x(6) - x(7)); (4/6.4)*(x(7) - x(8))];
f6 = @(t,x) [-k1(7)*x(1)*epo + 2*.10658*x(8); -k2(7)*(x(2))^2 + k1(7)*x(1)*epo; -.1066*x(3) + .5*k2(7)*((x(2))^2); .1066*x(3) - .10658*x(8); (4/6.4)*(x(3) - x(5)); (4/6.4)*(x(5) - x(6)); (4/6.4)*(x(6) - x(7)); (4/6.4)*(x(7) - x(8))];
f7 = @(t,x) [-k1(8)*x(1)*epo + 2*.10658*x(8); -k2(8)*(x(2))^2 + k1(8)*x(1)*epo; -.1066*x(3) + .5*k2(8)*((x(2))^2); .1066*x(3) - .10658*x(8); (4/6.4)*(x(3) - x(5)); (4/6.4)*(x(5) - x(6)); (4/6.4)*(x(6) - x(7)); (4/6.4)*(x(7) - x(8))];
f8 = @(t,x) [-k1(9)*x(1)*epo + 2*.10658*x(8); -k2(9)*(x(2))^2 + k1(9)*x(1)*epo; -.1066*x(3) + .5*k2(9)*((x(2))^2); .1066*x(3) - .10658*x(8); (4/6.4)*(x(3) - x(5)); (4/6.4)*(x(5) - x(6)); (4/6.4)*(x(6) - x(7)); (4/6.4)*(x(7) - x(8))];


x0 = [1 0 0 0 0 0 0 0];
options = odeset('NonNegative',1);
[t, xno0] = ode15s(f0,[0:.1:5],x0,options);
[t, xno1] = ode15s(f1,[0:.1:5],x0,options);
[t, xno2] = ode15s(f2,[0:.1:5],x0,options);
[t, xno3] = ode15s(f3,[0:.1:5],x0,options);
[t, xno4] = ode15s(f4,[0:.1:5],x0,options);
[t, xno5] = ode15s(f5,[0:.1:5],x0,options);
[t, xno6] = ode15s(f6,[0:.1:5],x0,options);
[t, xno7] = ode15s(f7,[0:.1:5],x0,options);
[t, xno8] = ode15s(f8,[0:.1:5],x0,options);
xno = [xno0 xno1 xno2 xno3 xno4 xno5 xno6 xno7 xno8];

x1exval = pi*[xno(:,1) xno(:,9) xno(:,17) xno(:,25) xno(:,33) xno(:,41) xno(:,49) xno(:,57) xno(:,65)]';
x2exval = pi*[xno(:,2) xno(:,10) xno(:,18) xno(:,26) xno(:,34) xno(:,42) xno(:,50) xno(:,58) xno(:,66)]';
x3exval = pi*[xno(:,3) xno(:,11) xno(:,19) xno(:,27) xno(:,35) xno(:,43) xno(:,51) xno(:,59) xno(:,67)]';
x4exval = pi*[xno(:,4) xno(:,12) xno(:,20) xno(:,28) xno(:,36) xno(:,44) xno(:,52) xno(:,60) xno(:,68)]';

x1exnew = [xno(:,1) xno(:,9) xno(:,17) xno(:,25) xno(:,33) xno(:,41) xno(:,49) xno(:,57) xno(:,65)]';
x2exnew = [xno(:,2) xno(:,10) xno(:,18) xno(:,26) xno(:,34) xno(:,42) xno(:,50) xno(:,58) xno(:,66)]';
x3exnew = [xno(:,3) xno(:,11) xno(:,19) xno(:,27) xno(:,35) xno(:,43) xno(:,51) xno(:,59) xno(:,67)]';
x4exnew = [xno(:,4) xno(:,12) xno(:,20) xno(:,28) xno(:,36) xno(:,44) xno(:,52) xno(:,60) xno(:,68)]';



save('x1exnew','x1exnew')
save('x2exnew','x2exnew')
save('x3exnew','x3exnew')
save('x4exnew','x4exnew')
save('x1eval');
save('x2eval');
save('x3eval');
save('x4eval');

%to plot mean, do plot(t,x1mean(1,:)), etc etc
figure
plot(t,x1exval(1,:)); xlabel('Time (min)'); ylabel('Expected Value');
hold on; plot(t,x2exval(1,:)); plot(t,x3exval(1,:)); plot(t,x4exval(1,:))
legend('X_{1}','X_{2}','X_{3}','X_{4}')
set(gca,'FontSize',15)

for i = 1:1:6
    var1(i,:) = [((x1exval(i,:)).^2).*phiprod(i)];
    var2(i,:) = [((x2exval(i,:)).^2).*phiprod(i)];
    var3(i,:) = [((x3exval(i,:)).^2).*phiprod(i)];
    var4(i,:) = [((x4exval(i,:)).^2).*phiprod(i)];
end

var1 = sum(var1) - (x1exval(1,:)).^2;
var2 = sum(var2) - (x2exval(1,:)).^2;
var3 = sum(var3) - (x3exval(1,:)).^2;
var4 = sum(var4) - (x4exval(1,:)).^2;

figure 
plot(t,var1); xlabel('Time (min)'); ylabel('Variance');
hold on;
plot(t,var2); 
plot(t,var3); 
plot(t,var4); 
legend('X_{1}','X_{2}','X_{3}','X_{4}')
set(gca,'FontSize',15)