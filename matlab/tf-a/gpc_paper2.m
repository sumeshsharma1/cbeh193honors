function dx = gpc_paper2(t,x0,flag,omega)

load('D1');
load('D2');
load('D3');
load('E');
load('a');
load('C');
load('omega');

x1_coeff = x0(1:3)';
q1_coeff = x0(4:6)';
q2_coeff = x0(7:9)';
q3_coeff = x0(10:12)';
q4_coeff = x0(13:15)';
q5_coeff = x0(16:18)';
q6_coeff = x0(19:21)';
q7_coeff = x0(22:24)';
q8_coeff = x0(25:27)';

dx1 = ((20*kron(q8_coeff,q8_coeff)*D1*E)/(x1_coeff(1).^2 + 10)) - (0.1*x1_coeff) + 0.01;
dq1 = 8*(kron(a,x1_coeff) - kron(a,q1_coeff))*D1*E;
dq2 = 8*(kron(a,q1_coeff) - kron(a,q2_coeff))*D1*E;
dq3 = 8*(kron(a,q2_coeff) - kron(a,q3_coeff))*D1*E;
dq4 = 8*(kron(a,q3_coeff) - kron(a,q4_coeff))*D1*E;
dq5 = 8*(kron(a,q4_coeff) - kron(a,q5_coeff))*D1*E;
dq6 = 8*(kron(a,q5_coeff) - kron(a,q6_coeff))*D1*E;
dq7 = 8*(kron(a,q6_coeff) - kron(a,q7_coeff))*D1*E;
dq8 = 8*(kron(a,q7_coeff) - kron(a,q8_coeff))*D1*E;

dx = [dx1 dq1 dq2 dq3 dq4 dq5 dq6 dq7 dq8]';


