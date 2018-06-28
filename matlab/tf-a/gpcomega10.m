n_xi = 2;
m = 1;
P = factorial(n_xi+m)/(factorial(n_xi)*factorial(m));
omega = 1;
for k=1:P
x1(k) = sym(sprintf('x1%d(t)', k));
end
x = x1;
for j = 1:n_xi
    xi(j) = sym(sprintf('xi%d', j));
end
phi1 = 1; 
phi2 = xi(1);
%phi3 = xi(2);
phi3 = (xi(1))^2 - 1;
phi5 = (xi(1))*(xi(2));
phi6 = (xi(2))^2 - 1;
C = [phi1 phi2 phi3];
for k = 1:P
    phiprod(k) = inner_product_multi_d('hermite',C(k),C(k),xi);
end

D1 = inner_product_multi_d('hermite',kron(C,C)',C,xi);
D2 = inner_product_multi_d('hermite',kron(kron(C,C),C)',C,xi);
D3 = inner_product_multi_d('hermite',kron(kron(kron(C,C),C),C)',C,xi);
E = diag(1./phiprod);
for k = 1:P
     %a(k) = (inner_product_multi_d('hermite',((1.176e-4)*xi(1) + (1/120)),C(k),xi))/(phiprod(k)); %Gaussian
     a(k) = (inner_product_multi_d('legendre',(1/(120 + xi(1)*10)),C(k),xi))/(phiprod(k)); %Legendre
end

save('D1','D1');
save('D2','D2');
save('D3','D3');
save('E','E');
save('a','a');
save('C','C');
save('omega','omega');
save('phiprod','phiprod');


