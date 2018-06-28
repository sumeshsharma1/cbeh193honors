load('vartotal1'); load('vartotal2'); load('vartotal3'); load('vartotal4');
load('x1exnew'); load('x2exnew'); load('x3exnew'); load('x4exnew'); 
load('pi'); load('phiprod')

for i = 1:1:9 %N is 9, L is 6
    varsum1(i,:) = pi(1,i).*vartotal1(i,:); %for each state
    varsum2(i,:) = pi(1,i).*vartotal2(i,:);
    varsum3(i,:) = pi(1,i).*vartotal3(i,:);
    varsum4(i,:) = pi(1,i).*vartotal4(i,:);
end
tmp1real = zeros(1,51);
tmp2real = zeros(1,51);
tmp3real = zeros(1,51);
tmp4real = zeros(1,51);
for i = 1:1:9
    for j = 1:1:9
        for k = 1:1:6
          tmp1 = pi(k,i)*pi(k,j)*phiprod(k)*x1exnew(i,:).*x1exnew(j,:);
          tmp1real = tmp1real+tmp1;
          tmp2 = pi(k,i)*pi(k,j)*phiprod(k)*x2exnew(i,:).*x2exnew(j,:);
          tmp2real = tmp2real+tmp2;
          tmp3 = pi(k,i)*pi(k,j)*phiprod(k)*x3exnew(i,:).*x3exnew(j,:);
          tmp3real = tmp3real+tmp3;
          tmp4 = pi(k,i)*pi(k,j)*phiprod(k)*x4exnew(i,:).*x4exnew(j,:);
          tmp4real = tmp4real+tmp4;
        end   
    end
end

totvar1 = sum(varsum1,1) + sum(tmp1real,1) - (pi(1,:)*x1exnew).^2;
totvar2 = sum(varsum2,1) + sum(tmp2real,1) - (pi(1,:)*x2exnew).^2;
totvar3 = sum(varsum3,1) + sum(tmp3real,1) - (pi(1,:)*x3exnew).^2;
totvar4 = sum(varsum4,1) + sum(tmp4real,1) - (pi(1,:)*x4exnew).^2;

t = [0:.1:5];

save('totvar1','totvar1')
save('totvar2','totvar2')
save('totvar3','totvar3')
save('totvar4','totvar4')

figure
hold on 
plot(t,totvar1), plot(t,totvar2),plot(t,totvar3), plot(t,totvar4)
xlabel('Time')
ylabel('Variance')
legend('X1','X2','X3','X4')
title('Jak/STAT Total Variance')