function dx = zebra_paper(t,x0)

load('D1');
load('E');
load('Tp');
load('C');
load('Tm');
load('phinew');
load('phiprod');

p1_coeff = x0(1:6)';
m1_coeff = x0(7:12)';
qp1_coeff = x0(13:18)';
qp2_coeff = x0(19:24)';
qp3_coeff = x0(25:30)';
qp4_coeff = x0(31:36)';
qm1_coeff = x0(37:42)';
qm2_coeff = x0(43:48)';
qm3_coeff = x0(49:54)';
qm4_coeff = x0(55:60)';

dp1 = 4.5*qp4_coeff - 0.23*p1_coeff;
dm1 = (33/(1 + (qm4_coeff(1).^2)/100))*(phinew./phiprod) - 0.23*m1_coeff;

dqp1 = 4*(kron(Tp,m1_coeff) - kron(Tp,qp1_coeff))*D1*E;
dqp2 = 4*(kron(Tp,qp1_coeff) - kron(Tp,qp2_coeff))*D1*E;
dqp3 = 4*(kron(Tp,qp2_coeff) - kron(Tp,qp3_coeff))*D1*E;
dqp4 = 4*(kron(Tp,qp3_coeff) - kron(Tp,qp4_coeff))*D1*E;
dqm1 = 4*(kron(Tm,p1_coeff) - kron(Tm,qm1_coeff))*D1*E;
dqm2 = 4*(kron(Tm,qm1_coeff) - kron(Tm,qm2_coeff))*D1*E;
dqm3 = 4*(kron(Tm,qm2_coeff) - kron(Tm,qm3_coeff))*D1*E;
dqm4 = 4*(kron(Tm,qm3_coeff) - kron(Tm,qm4_coeff))*D1*E;

dx = [dp1 dm1 dqp1 dqp2 dqp3 dqp4 dqm1 dqm2 dqm3 dqm4]';

