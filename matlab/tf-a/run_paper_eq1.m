load('phiprod')
omega = 1;
t1 = [0:1:200];
t2 = [200:1:1000];
%R = normrnd(120,(10/6)^2,[1,length(t)]);
x0 = zeros(1,27);
options = odeset('NonNegative',1);
[t1, xvalues] = ode15s('gpc_paper',t1,x0,options,omega);

x1 = xvalues(length(t1),:);
[t2, xvalues2] = ode15s('gpc_paper2',t2,x1,options,omega);

varphi = [];
for i = 1:1:3
    varphi(i) = (phiprod(i)).^2;
    xvar(:,i) = ((xvalues(:,i)).^2)*varphi(i);
    xvar2(:,i) = ((xvalues2(:,i)).^2)*varphi(i);
end

xvar = sum(xvar,2) - ((xvalues(:,1)).^2);
xvar2 = sum(xvar2,2) - ((xvalues2(:,1)).^2);

    

