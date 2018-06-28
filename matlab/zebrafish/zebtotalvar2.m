load('vartotnew1'); load('vartotnew2'); load('vartotnew3'); load('vartotnew4')
load('vartotnew5'); load('vartotnew6'); load('vartotnew7'); load('vartotnew8')
load('vartotnew9'); load('x1znew'); load('x2znew'); load('pi'), load('phiprod')

varmat1 = [vartotnew1;vartotnew2;vartotnew3;vartotnew4;vartotnew5;vartotnew6;vartotnew7;vartotnew8;vartotnew9];

for i = 1:1:9 %N is 9, L is 6
    varsum1(i,:) = pi(1,i).*varmat1(2*i - 1,:); %State 1
    varsum2(i,:) = pi(1,i).*varmat1(2*i,:); %State 2
end
tmp1real = zeros(1,100001);
tmp2real = zeros(1,100001);
for i = 1:1:9
    for j = 1:1:9
        for k = 1:1:6
           tmp1 = pi(k,i)*pi(k,j)*phiprod(k)*x1znew(i,:).*x1znew(j,:);
           tmp1real = tmp1real+tmp1;
           tmp2 = pi(k,i)*pi(k,j)*phiprod(k)*x2znew(i,:).*x2znew(j,:);
           tmp2real = tmp2real+tmp2;
        end
    end
end

% Total variance for each state
totvar1 = sum(varsum1,1) + sum(tmp1real,1) - (pi(1,:)*x1znew).^2;
totvar2 = sum(varsum2,1) + sum(tmp2real,1) - (pi(1,:)*x2znew).^2;

t = [0:0.01:1000];

figure
hold on 
plot(t,totvar1), plot(t,totvar2), ylim([0 inf])
xlabel('Time')
ylabel('Variance')
legend('P(t)','M(t)')
title('Zebrafish Variance for EKF+gPC')
        

