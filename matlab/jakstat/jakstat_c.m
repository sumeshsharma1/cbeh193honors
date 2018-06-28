function dx = jakstat_c(t,x0,flag,epo)
load('D1');
load('D2');
load('E');
load('k1');
load('k2');
load('C');
k3 = .1066;
k4 = .10658;

x1_coeff = x0(1:6)';
x2_coeff = x0(7:12)';
x3_coeff = x0(13:18)';
x4_coeff = x0(19:24)';
q1_coeff = x0(25:30)';
q2_coeff = x0(31:36)';
q3_coeff = x0(37:42)';
q4_coeff = x0(43:48)';


dx1 = -((kron(k1,x1_coeff))*D1*E*epo + 2*k4*q4_coeff);
dx2 = -(kron(kron(k2,x2_coeff), x2_coeff)*D2*E) + kron(k1,x1_coeff)*D1*E.*epo;
dx3 = -k3*x3_coeff + .5*(kron(kron(k2,x2_coeff), x2_coeff)*D2*E);
dx4 = k3*x3_coeff - k4*q4_coeff;
dq1 = (4/6.4)*(x3_coeff - q1_coeff);
dq2 = (4/6.4)*(q1_coeff - q2_coeff);
dq3 = (4/6.4)*(q2_coeff - q3_coeff);
dq4 = (4/6.4)*(q3_coeff - q4_coeff);

dx = [dx1 dx2 dx3 dx4 dq1 dq2 dq3 dq4]';

