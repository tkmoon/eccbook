function z = getinfs(data)

bins = -100:0.1:100;
for n=1:size(data,1)
    [N(n,:),X(n,:)] = hist(log((data(n,:)+realmin) ./ (1-data(n,:)+realmin)),bins);
    [I(n)] = getinf(N(n,:),X(n,:));
end

% figure(1)
% plot(0:n-1,I,'x-');
% xlabel('Iteration Number');
% ylabel('Mutual Information');
% title('Mutual Information path');
% 
% figure(2)
% plot(X(1,:),Ib(1,:),'r');
% hold on;
% plot(X(end,:),Ib(end,:),'b');
% legend('Received Signal','Final Iteration');
% xlabel('Log likelihood ratio');
% ylabel('Frequency');
% title('Initial and Final LL Distributions');

z = I;