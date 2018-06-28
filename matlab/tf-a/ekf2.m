%% EKF system: Change Tau for different results, reference ekf.m
t1 = [0:1:200];
load('xp1.mat')
load('Ppic')
initial_time = 201;
final_time = 1000;
tstep = 1;
Nx = 5;
Ny = 1;
N_sim = (final_time - initial_time)/tstep +1;

t = initial_time:tstep:final_time;
%initialise all vectors as a function of NaN
%Vectors of states and outputs:
x = zeros(Nx,N_sim); % X S P
y = zeros(Ny,N_sim);
xp = zeros(Nx,N_sim);
%Noises:
w = NaN(Nx,Nx,N_sim); %Process noise
v = NaN(Ny,Ny,N_sim); %Measurement noise

% Matrices for noise:
%Q = eye(Nx)*0.01;
%R = eye(Ny)*0.01;

%Initialise parameters:
%Km = 1.2;
%mu_max = 0.48;
%alpha = 2.2;
%beta = 0.2;
%Sf = 20;
%Yxs = 0.4;

omega = 1;
kd1 = 0.1;
R_bas = .01*omega;
Kd = 10*omega^2;
kf = 20*omega;   
tau = tau; %tau varied between .0085, 1/120, and .0081



% dX = (-D*X + mu*X); % continuous x(1,k)
% dS = (D*(Sf - S) - (1/Yxs)*mu*X); 
% dP = (-D*P + ((alpha*mu)+beta)*X);
    %No inhibition:
%     mu = (mu_max*S)/(Km+S);

S = zeros(Ny,Ny,N_sim);
invS = zeros(Ny,Ny,N_sim);
K = zeros(Nx,Ny,N_sim);
P = zeros(Nx,Nx,N_sim);
Pp = zeros(Nx,Nx,N_sim);
ytilde = zeros(Ny,N_sim);
mu = zeros(1,N_sim);

% initialise the whole thing with time step = 1;
x(1:5,1) = [xp1(1,length(xp1(1,:))) xp1(2,length(xp1(2,:))) xp1(3,length(xp1(3,:))) xp1(4,length(xp1(4,:))) xp1(5,length(xp1(5,:)))]';
%mu_max = 0.48;
%x(4,1) = 0.48;
dFdx(:,:,1) = zeros(5,5);
dHdx(:,:,1) = [1 0 0 0 0]';
%mu(1) = mu_max*x(2,1)/(Km+x(2,1));
%D = zeros(1,N_sim);
%D(1,1) = 0.15;
xp(1:5,1) = x(:,1);
z(:,1) = xp(1,1);

%% Inputs:
% D(1,1) = 
% Sf(1,1) =

for k = 2:N_sim %either this or from 1:N_sim-1 taking into consideration the first iteration \
    k;
%just xp(1,k)
    
    xp(1,k) = xp(1,k-1)+ (((kf*((x(5,k-1))^2))/((x(5,k-1))^2 + Kd)) - kd1*x(1,k-1) + R_bas)*tstep;
    xp(2,k) = xp(2,k-1) + tstep*((4*tau)*(x(1,k-1) - x(2,k-1)));
    xp(3,k) = xp(3,k-1) + tstep*((4*tau)*(x(2,k-1) - x(3,k-1)));
    xp(4,k) = xp(4,k-1) + tstep*((4*tau)*(x(3,k-1) - x(4,k-1)));
    xp(5,k) = xp(5,k-1) + tstep*((4*tau)*(x(4,k-1) - x(5,k-1)));
    %xp(1,k)  = x(1,k-1) + tstep*(-D(k-1)*x(1,k-1) + mu(k)*x(1,k-1)); %discrete x 

    
    % Jacobian analytical expressions: - there needs to be an extra dimension
    % because of the uncertain parameter:
    
% Use delay equations for time delay term to get 5x5 matrix
% (f1,f2,f3,f4,f5), take derivatives wrt x1,x2,x3,x4,x5
    dFdx(:,:,k) = tstep*[[-kd1,0,0,0,((2*kf*Kd*x(5,k-1))/((Kd + ((x(5,k-1))^2))^2))];
        [4*tau, -4*tau, 0,0,0]; [0,4*tau,-4*tau,0,0]; [0,0,4*tau,-4*tau,0]; [0,0,0,4*tau,-4*tau]];
    dHdx(:,:,k) = [1 0 0 0 0];
    Q = eye(Nx)*sqrt((kf*(x(5,k-1))^2)/((x(5,k-1))^2 + Kd) + kd1*x(1,k-1)); %multiply by noise coeff
    R = 0;
    Pp(:,:,k) = dFdx(:,:,k)*P(:,:,k-1)*dFdx(:,:,k)' + Q;
    %S(:,:,k) = dHdx(:,:,k)*Pp(:,:,k)*dHdx(:,:,k)' + R;
    %invS(:,:,k) = inv(S(:,:,k));
    
    %K(:,:,k) = Pp(:,:,k)*dHdx(:,:,k)'*(invS(:,:,k)); 
    
    %Measurement residual:
    z(:,k) = [xp(1,k)];
    %ytilde(:,k) = z(:,k) - [xp(1,k)]; 
    %Update the state estimate as well as the covariance estimate:
    x(:,k) = xp(:,k);% + K(:,:,k)*ytilde(:,k);
    %P(:,:,k) = (eye(Nx) - K(:,:,k)*dHdx(:,:,k))*Pp(:,:,k);
    P(:,:,k) = Pp(:,:,k);
    A(:,:,k) = (dFdx(:,:,k)*tstep + eye(Nx)); % [[-D(k-1) + mu(k),0,0,0];[-mu(k)/Yxs , D(k-1)*Sf/x(2,k-1) - D(k-1),0,0];[alpha*mu(k)+beta, 0 ,-D(k-1),0];[0, 0 , 0 , 1]];
    %C(:,:,k) = [1 0 0 0 ; 0 1 0 0 ];
    %D(k) = 0.01+D(k-1);
    
 %plot(1,1) of Pp matrix for variance   
end

figure 
subplot(1,2,1),
plot(t,xp(1,:)); hold on; plot(t1,xp1(1,:))
xlabel('Time'); ylabel('Expected Value'); title(['Expected Value for Tau = ',num2str(tau)])

subplot(1,2,2),
plot(t,squeeze(Pp(1,1,:))'); hold on; plot(t1,Ppic)
xlabel('Time'); ylabel('Variance'); title(['Variance for Tau = ',num2str(tau)])

var0085 = [Ppic, squeeze(Pp(1,1,:))'];

save('var0085','var0085')
