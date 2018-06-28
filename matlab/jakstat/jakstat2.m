n_xi = 2;
m = 2;
P = factorial(n_xi+m)/(factorial(n_xi)*factorial(m));
k3 = .1066;
k4 = .10658;
for k=1:P
x1(k) = sym(sprintf('x1%d(t)', k));
x2(k) = sym(sprintf('x2%d(t)', k));
x3(k) = sym(sprintf('x3%d(t)', k));
x4(k) = sym(sprintf('x4%d(t)', k));
end
x = [x1 x2 x3 x4];
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
end
D1 = inner_product_multi_d('hermite',kron(C,C)',C,xi);
D2 = inner_product_multi_d('hermite',kron(kron(C,C),C)',C,xi);
E = diag(1./phiprod);
for k = 1:P
    k1(k) = (inner_product_multi_d('hermite',(.1*xi(1) + .21),C(k),xi))/(phiprod(k));
    k2(k) = (inner_product_multi_d('hermite',(.1*xi(2) + 2.46),C(k),xi))/(phiprod(k));
    x3new(k) = (x(12+k))*phiprod(k);
end

save('k1','k1')
save('k2','k2')
