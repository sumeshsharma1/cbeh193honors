%% EKF system

final_time = 1000;
tstep = 0.1;
Nx = 10;
Ny = 1;
N_sim = final_time/tstep +1;

t = 0:tstep:final_time;
%initialise all vectors as a function of NaN
%Vectors of states and outputs:
x = zeros(Nx,N_sim); % X S P
y = zeros(Ny,N_sim);
xp = zeros(Nx,N_sim);
%Noises:
w = NaN(Nx,Nx,N_sim); %Process noise
v = NaN(Ny,Ny,N_sim); %Measurement noise


%Initialise parameters:

axi1 = [-sqrt(3) -sqrt(3) -sqrt(3) 0 0 0 sqrt(3) sqrt(3) sqrt(3)]; 
axi2 = [-sqrt(3) 0 sqrt(3) -sqrt(3) 0 sqrt(3) -sqrt(3) 0 sqrt(3)];

Tm = 12;
Tp = 2.8;
a = .225;
m = 33;
b = 0.23;
c = 0.23;
p0 = 10;
vartotal1 = [];
vartotal2 = [];

%for i = 1:1:9
tau1 = .28*axi1(5) + (2.8); %Tp
tau2 = 1.2*axi2(5) + (12); %Tm



S = zeros(Ny,Ny,N_sim);
invS = zeros(Ny,Ny,N_sim);
K = zeros(Nx,Ny,N_sim);
P = zeros(Nx,Nx,N_sim);
Pp = zeros(Nx,Nx,N_sim);
ytilde = zeros(Ny,N_sim);
mu = zeros(1,N_sim);

% initialise the whole thing with time step = 1;
%p(t) is first, then m(t), then p(t) time delay, then m(t) time delay
x(1:10,1) = [10 1 0 0 0 0 0 0 0 0]';
dFdx(:,:,1) = zeros(10,10);
dHdx(:,:,1) = [1 1 0 0 0 0 0 0 0 0]';
xp(1:10,1) = x(:,1);
z(:,1) = xp(1,1);

%% Inputs:

for k = 2:N_sim %either this or from 1:N_sim-1 taking into consideration the first iteration \
    %replace k with m, do jacobian with tstep
%just xp(1,k)
    
    xp(1,k) = xp(1,k-1) + (a*x(10,k-1) - b*x(1,k-1))*tstep; %p(t)
    xp(2,k) = xp(2,k-1) + tstep*((m/(1 + ((x(6,k-1))^2)/(p0^2))) - c*x(2,k-1)); %m(t)
    xp(3,k) = xp(3,k-1) + tstep*((4/tau2)*(x(1,k-1) - x(3,k-1))); %p(t) time delay
    xp(4,k) = xp(4,k-1) + tstep*((4/tau2)*(x(3,k-1) - x(4,k-1)));
    xp(5,k) = xp(5,k-1) + tstep*((4/tau2)*(x(4,k-1) - x(5,k-1)));
    xp(6,k) = xp(6,k-1) + tstep*((4/tau2)*(x(5,k-1) - x(6,k-1)));
    xp(7,k) = xp(7,k-1) + tstep*((4/tau1)*(x(2,k-1) - x(7,k-1))); %m(t) time delay
    xp(8,k) = xp(8,k-1) + tstep*((4/tau1)*(x(7,k-1) - x(8,k-1)));
    xp(9,k) = xp(9,k-1) + tstep*((4/tau1)*(x(8,k-1) - x(9,k-1)));
    xp(10,k) = xp(10,k-1) + tstep*((4/tau1)*(x(9,k-1) - x(10,k-1)));
    

    
    % Jacobian analytical expressions: - there needs to be an extra dimension
    % because of the uncertain parameter:
    
% Use delay equations for time delay term to get 10x10 matrix
    dFdx(:,:,k) = [[-b,0,0,0,0,0,0,0,0,a];
        [0,-c,0,0,0,-2*m*(p0^2)*x(6,k-1)/((p0^2 + (x(6,k-1))^2)^2),0,0,0,0]; 
        [4/tau2,0,-4/tau2,0,0,0,0,0,0,0]; [0,0,4/tau2,-4/tau2,0,0,0,0,0,0];
        [0,0,0,4/tau2,-4/tau2,0,0,0,0,0]; [0,0,0,0,4/tau2,-4/tau2,0,0,0,0];
        [0,4/tau1,0,0,0,0,-4/tau1,0,0,0]; [0,0,0,0,0,0,4/tau1,-4/tau1,0,0];
        [0,0,0,0,0,0,0,4/tau1,-4/tau1,0]; [0,0,0,0,0,0,0,0,4/tau1,-4/tau1]];
    dHdx(:,:,k) = [1 1 0 0 0 0 0 0 0 0];
    Q = zeros(Nx); 
    Q(1,1) = (sqrt(abs(a*x(10,k-1))) - sqrt(abs(b*x(1,k-1))))^2;
    %Q(1,1) = abs(a*x(10,k-1) - (b*x(1,k-1))); %these two are the real 1s
    Q(2,2) = (sqrt((m/(1 + ((x(6,k-1))^2)/(p0^2)))) - sqrt(abs(c*x(2,k-1))))^2; %multiply by noise coeff
    %Q(2,2) = abs((m/(1 + ((x(6,k-1))^2)/(p0^2))) - c*x(2,k-1));
%     for j = 1:1:2 %Figure 5
%         Q(j,j) = tstep*(.005)*xp(j,k-1);
%     end
%     Q(1,1) = .01;
%     Q(2,2) = .01;
    R = 0;
    %Pp(:,:,k) = dFdx(:,:,k)*P(:,:,k-1)*dFdx(:,:,k)' + tstep*Q;
    Pp(:,:,k) = Pp(:,:,k-1) + tstep*(dFdx(:,:,k)*P(:,:,k-1) + P(:,:,k-1)*dFdx(:,:,k)' + Q);
    
    %Measurement residual:
    z(:,k) = [xp(1,k)];
    %Update the state estimate as well as the covariance estimate:
    x(:,k) = xp(:,k);
    P(:,:,k) = Pp(:,:,k);
    A(:,:,k) = (dFdx(:,:,k)*tstep + eye(Nx)); % [[-D(k-1) + mu(k),0,0,0];[-mu(k)/Yxs , D(k-1)*Sf/x(2,k-1) - D(k-1),0,0];[alpha*mu(k)+beta, 0 ,-D(k-1),0];[0, 0 , 0 , 1]];
   
    
 %plot(1,1) of Pp matrix for variance   
end
var1 = squeeze(Pp(1,1,:));
var2 = squeeze(Pp(2,2,:));
% vartotal1(i,:) = [var1]';
% vartotal2(i,:) = [var2]';
% end 
% 
% save('vartotal1','vartotal1')
% save('vartotal2','vartotal2')


figure
% subplot(1,2,1),
% plot(t,xp(1,:)); hold on; plot(t,xp(2,:));
% xlabel('Time'), ylabel('Expected Value');
% title(['Expected Value for Tm = ',num2str(tau2,3), ' and Tp = ',num2str(tau1,3)])
% legend('p(t)','m(t)')
% 
% subplot(1,2,2),
plot(t,var1); hold on; plot(t,var2);
xlabel('Time'), ylabel('Variance');
legend('p(t)','m(t)')
set(gca,'fontsize',15)



