% this file studies the second eigenvalue of an Erdos-Renyi random graph

clear all; close all; %rng(2017);

%% parameter

n = 100;
n = 200;
n = 500;

%n = 1000;
%n = 5000;

%

%a = 4;
%p = a*log(n)/n

p = 0.15;

%% simulate random graph

W = double(triu( rand(n,n) < p, 1));

W= W + W';

%
figure(1),clf;
%% the degree 
dW = sum(W,2);

imagesc(W);
title('W')


figure(2),clf; hold on;
plot(dW,'x-');
plot( p*n*ones(n,1),'--');
grid on;
xlabel('node index')
legend('degree of node', 'expectation')
title('degree')


%% eig

[v,d]=eig(W,diag(dW));
[lambda, tmp] = sort(diag(d),'descend');
psi = v(:,tmp);

Lambda_L = 1-lambda;

%
figure(3),clf;
subplot(2,1,1);
plot(Lambda_L,'x-');
axis([1,n,0,2])
grid on; title('Lambda of L_{sym}')
subplot(2,1,2);
hist(Lambda_L, 100);
grid on; title('histogram')

%%
sig= sqrt((1-p)/p); %due to normalize by degree

C = 2*sig


figure(4),clf; hold on;
%hist((1-Lambda_L(2:n))*sqrt(n),100);
[n1,x1]=hist((1-Lambda_L(2:n))*sqrt(n),100);
bar(x1,n1,'FaceColor','w','EdgeColor','b');
plot(-C,0,'xr', 'MarkerSize',10);
plot(C,0,'xr', 'MarkerSize',10);
grid on; 
title('histogram of 1-Lambda(2:n)*sqrt(n)')



