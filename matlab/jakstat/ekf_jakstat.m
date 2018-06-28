%% EKF system: NOT FULLY DONE FOR SYSTEM YET 

final_time = 5;
tstep = .1;
Nx = 8;
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

% Matrices for noise:
%Q = eye(Nx)*0.01;
%R = eye(Ny)*0.01;
xi1 = [-sqrt(3) -sqrt(3) -sqrt(3) 0 0 0 sqrt(3) sqrt(3) sqrt(3)];
xi2 = [-sqrt(3) 0 sqrt(3) -sqrt(3) 0 sqrt(3) -sqrt(3) 0 sqrt(3)];
%Initialise parameters:
for i = 1:1:9
k1 = xi1(i)*.0021 + .021;
k2 = xi2(i)*.246 + 2.46;
k3 = .1066; %min^-1
k4 = .10658; %min^-1
tau = 6.4; %min
epo = 100;

S = zeros(Ny,Ny,N_sim);
invS = zeros(Ny,Ny,N_sim);
K = zeros(Nx,Ny,N_sim);
P = zeros(Nx,Nx,N_sim);
Pp = zeros(Nx,Nx,N_sim);
ytilde = zeros(Ny,N_sim);
mu = zeros(1,N_sim);

% initialise the whole thing with first concentration 1, rest 0 
x(1:8,1) = [1 0 0 0 0 0 0 0]';
dFdx(:,:,1) = zeros(8,8);
dHdx(:,:,1) = [1 1 1 1 0 0 0 0]';
xp(1:8,1) = x(:,1);
z(:,1) = xp(1,1);

%% Inputs:
% D(1,1) = 
% Sf(1,1) =

for k = 2:N_sim %either this or from 1:N_sim-1 taking into consideration the first iteration \
    
%just xp(1,k)
    
    xp(1,k) = xp(1,k-1) + tstep*(-k1*x(1,k-1)*epo + 2*k4*x(8,k-1));
    xp(2,k) = xp(2,k-1) + tstep*(-k2*((x(2,k-1))^2) + k1*x(1,k-1)*epo);
    xp(3,k) = xp(3,k-1) + tstep*(-k3*x(3,k-1) + .5*k2*((x(2,k-1))^2));
    xp(4,k) = xp(4,k-1) + tstep*(k3*x(3,k-1) - k4*x(8,k-1));
    xp(5,k) = xp(5,k-1) + tstep*((4/tau)*(x(3,k-1) - x(5,k-1)));
    xp(6,k) = xp(6,k-1) + tstep*((4/tau)*(x(5,k-1) - x(6,k-1)));
    xp(7,k) = xp(7,k-1) + tstep*((4/tau)*(x(6,k-1) - x(7,k-1))); 
    xp(8,k) = xp(8,k-1) + tstep*((4/tau)*(x(7,k-1) - x(8,k-1)));
  
    

    
    % Jacobian analytical expressions: - there needs to be an extra dimension
    % because of the uncertain parameter:
    
% Use delay equations for time delay term to get 8X8 matrix
    dFdx(:,:,k) = [[-k1*epo,0,0,0,0,0,0,2*k4];
        [k1*epo,-2*k2*x(2,k-1),0,0,0,0,0,0]; 
        [0,k2*x(2,k-1),-k3,0,0,0,0,0]; [0,0,k3,0,0,0,0,-k4];
        [0,0,4/tau,0,-4/tau,0,0,0]; [0,0,0,0,4/tau,-4/tau,0,0];
        [0,0,0,0,0,4/tau,-4/tau,0]; [0,0,0,0,0,0,4/tau,-4/tau]];
    dHdx(:,:,k) = [1 1 1 1 0 0 0 0];
    Q = zeros(Nx);
    if k == 2 
        Q = zeros(Nx);
    else
    c11 = var(xp(1,1:k-1));c12 = cov(xp(1,1:k-1),xp(2,1:k-1));
    c13 = cov(xp(1,1:k-1),xp(3,1:k-1)); c14 = cov(xp(1,1:k-1),xp(4,1:k-1));
    c22 = var(xp(2,1:k-1)); c23 = cov(xp(2,1:k-1),xp(3,1:k-1));
    c24 = cov(xp(2,1:k-1),xp(4,1:k-1)); c33 = var(xp(3,1:k-1));
    c34 = cov(xp(3,1:k-1),xp(4,1:k-1)); c44 = var(xp(4,1:k-1));
    %Figure 2
%     Q = [[c11,c12(1,2),c13(1,2),c14(1,2),0,0,0,0;]*xp(1,k-1)*.005
%             [c12(1,2),c22,c23(1,2),c24(1,2),0,0,0,0].*xp(2,k-1)*.005;
%             [c13(1,2),c23(1,2),c33,c34(1,2),0,0,0,0].*xp(3,k-1)*.005;
%             [c14(1,2),c24(1,2),c34(1,2),c44,0,0,0,0].*xp(4,k-1)*.005;
%             0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;
%             0,0,0,0,0,0,0,0;];
%       %Figure 3
%       Q = [[xp(1,k-1),sqrt(xp(2,k-1))*sqrt(xp(1,k-1)),sqrt(xp(3,k-1))*sqrt(xp(1,k-1)),sqrt(xp(4,k-1))*sqrt(xp(1,k-1)),0,0,0,0]*0.005;
%             [sqrt(xp(1,k-1))*sqrt(xp(2,k-1)),xp(2,k-1),sqrt(xp(3,k-1))*sqrt(xp(2,k-1)),sqrt(xp(4,k-1))*sqrt(xp(2,k-1)),0,0,0,0]*0.005;
%             [sqrt(xp(1,k-1))*sqrt(xp(3,k-1)),sqrt(xp(2,k-1))*sqrt(xp(3,k-1)),xp(3,k-1),sqrt(xp(4,k-1))*sqrt(xp(3,k-1)),0,0,0,0]*0.005;
%             [sqrt(xp(1,k-1))*sqrt(xp(4,k-1)),sqrt(xp(2,k-1))*sqrt(xp(4,k-1)),sqrt(xp(3,k-1))*sqrt(xp(4,k-1)),xp(4,k-1),0,0,0,0]*0.005;
%             0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;
%             0,0,0,0,0,0,0,0;];
%           
%       %Figure 4  
    Q = [[c11*xp(1,k-1),c12(1,2)*sqrt(xp(1,k-1)*xp(2,k-1)),c13(1,2)*sqrt(xp(1,k-1)*xp(3,k-1)),c14(1,2)*sqrt(xp(1,k-1)*xp(4,k-1)),0,0,0,0;]*.005
            [c12(1,2)*sqrt(xp(1,k-1)*xp(2,k-1)),c22*xp(2,k-1),c23(1,2)*sqrt(xp(2,k-1)*xp(3,k-1)),c24(1,2)*sqrt(xp(2,k-1)*xp(4,k-1)),0,0,0,0]*.005;
            [c13(1,2)*sqrt(xp(1,k-1)*xp(3,k-1)),c23(1,2)*sqrt(xp(2,k-1)*xp(3,k-1)),c33*xp(3,k-1),c34(1,2)*sqrt(xp(3,k-1)*xp(4,k-1)),0,0,0,0]*.005;
            [c14(1,2)*sqrt(xp(1,k-1)*xp(4,k-1)),c24(1,2)*sqrt(xp(2,k-1)*xp(4,k-1)),c34(1,2)*sqrt(xp(3,k-1)*xp(4,k-1)),c44*xp(4,k-1),0,0,0,0]*.005;
            0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0;];
     end
%     for j = 1:1:4 %Figure 5
%         Q(j,j) = tstep*(.005)*xp(j,k-1);
%     end
    R = 0;
    Pp(:,:,k) = tstep*(dFdx(:,:,k)*P(:,:,k-1)*dFdx(:,:,k)') + tstep*Q;
    %S(:,:,k) = dHdx(:,:,k)*Pp(:,:,k)*dHdx(:,:,k)' + R;
    %invS(:,:,k) = inv(S(:,:,k));
    
    %K(:,:,k) = Pp(:,:,k)*dHdx(:,:,k)'*(invS(:,:,k)); 
%     A(:,:,k) = (dFdx(:,:,k)*tstep + eye(Nx)); % [[-D(k-1) + mu(k),0,0,0];[-mu(k)/Yxs , D(k-1)*Sf/x(2,k-1) - D(k-1),0,0];[alpha*mu(k)+beta, 0 ,-D(k-1),0];[0, 0 , 0 , 1]];
%     Pp(:,:,k) = A(:,:,k)*P(:,:,k-1) + Q;
    %Measurement residual:
    z(:,k) = [xp(1,k)];
    %ytilde(:,k) = z(:,k) - [xp(1,k)]; 
    %Update the state estimate as well as the covariance estimate:
    x(:,k) = xp(:,k);% + K(:,:,k)*ytilde(:,k);
    %P(:,:,k) = (eye(Nx) - K(:,:,k)*dHdx(:,:,k))*Pp(:,:,k);
    P(:,:,k) = Pp(:,:,k);


 %plot(1,1) of Pp matrix for variance   
end

%figure 
% subplot(1,2,1),
% plot(t,xp(1,:)); hold on; plot(t,xp(2,:)); plot(t,xp(3,:)); plot(t,xp(4,:));
% xlabel('Time')
% ylabel('Expected Value')
% legend('X1','X2','X3','X4')
% title(['Expected Value for k1 = ',num2str(k1,3), ' and k2 = ',num2str(k2,3)])

var1 = squeeze(Pp(1,1,:));
var2 = squeeze(Pp(2,2,:));
var3 = squeeze(Pp(3,3,:));
var4 = squeeze(Pp(4,4,:));
vartotal1(i,:) = [var1]';
vartotal2(i,:) = [var2]';
vartotal3(i,:) = [var3]';
vartotal4(i,:) = [var4]';

end

save('vartotal1','vartotal1')
save('vartotal2','vartotal2')
save('vartotal3','vartotal3')
save('vartotal4','vartotal4')
figure
% subplot(1,2,2),
plot(t,var1); hold on; plot(t,var2); plot(t,var3); plot(t,var4);
xlabel('Time')
ylabel('Variance')
legend('X1','X2','X3','X4')
title(['Variance for k1 = ',num2str(k1,3), ' and k2 = ',num2str(k2,3)])

