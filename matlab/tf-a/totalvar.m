load('var0085') %load EKF variance for each sample
load('var0083')
load('var0081')
load('pi') %load nonintrusive pi matrix
load('a')
load('phiprod')
load('xnonew')
load('xeval')
load('mu')
varmat = [var0083;var0085;var0081];
for i = 1:1:3 %N is 3, L is 3
    varsum(i,:) = pi(1,i)*varmat(i,:);
end


tmpreal = zeros(1,1001);
for i = 1:1:3
    for j = 1:1:3
        for k = 1:1:3
           tmp = pi(k,i)*pi(k,j)*phiprod(k)*xnonew(i,:).*xnonew(j,:);
           tmpreal = tmpreal+tmp;
        end
    end
end
        
%for i = 1:1:4
%    itmp(i,:) = squeeze(jtmp(i,:,i));
%end

finalvar = sum(varsum,1) + sum(tmpreal,1) - mu.^2;
figure 
plot([0:1:1000],finalvar)
xlabel('Time')
ylabel('Variance')
title('Total Variance for 1 State Transcriptional Factor System')

    