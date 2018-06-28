n_xi = 2;
m = 2;
P = factorial(n_xi+m)/(factorial(n_xi)*factorial(m));
omega = 10;
for k=1:P
x1(k) = sym(sprintf('x1%d(t)', k));
end
x = x1;
for j = 1:n_xi
    xi(j) = sym(sprintf('xi%d', j));
end
phi1 = 1; 
phi2 = xi(1);
phi3 = xi(2);
phi4 = (xi(1))^2 - 1;
phi5 = (xi(1))*(xi(2));
phi6 = (xi(2))^2 - 1;
C = [phi1 phi2 phi3 phi4 phi5 phi6];
for k = 1:P
    phiprod(k) = inner_product_multi_d('hermite',C(k),C(k),xi);
    phinew(k) = inner_product_multi_d('hermite',C(k),1,xi);
end

D1 = inner_product_multi_d('hermite',kron(C,C)',C,xi);
E = diag(1./phiprod);
for k = 1:P
    %kft(k) = (inner_product_multi_d('hermite',(.01*xi(1) + 20*omega),C(k),xi))/(phiprod(k));%kf after t = 200
    Tp(k) = (inner_product_multi_d('hermite',((1.17e-3)*xi(1) + (1/12)),C(k),xi))/(phiprod(k));
    Tm(k) = (inner_product_multi_d('hermite',((4.8e-3)*xi(2) + (1/2.8)),C(k),xi))/(phiprod(k));
end

save('D1');
save('E');
save('Tp');
save('Tm');
save('C');
save('phiprod');
save('phinew');