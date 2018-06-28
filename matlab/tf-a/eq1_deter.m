function dx = eq1_deter(t,x0,flag,omega)

x1_coeff = x0(1:3)';
q1_coeff = x0(4:6)';

dx1 = (q1_coeff).^2 + 10*(omega^2);
dq1 = (1/120)*(x1_coeff-q1_coeff);

dx = [dx1 dq1]';