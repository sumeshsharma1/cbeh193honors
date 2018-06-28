for i = 1:1:9
    for j = 1:1:9
        tmp2 = zeros(1,501);
        for k = 1:6
           tmp = pi(k,i)*pi(k,j)*phiprod(k)*x2exnew(i,:).*x2exnew(j,:);
           tmp2 = tmp2+tmp;
        end
    end
end

plot([0:0.01:5],tmp2 + sum(vtmp2,1) - (pi(1,:)*x2exnew).^2)