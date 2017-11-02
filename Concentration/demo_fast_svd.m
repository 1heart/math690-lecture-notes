

clear all; close all; rng(2017)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% singluar values of gaussian matrix

n = 1000;

%k = 20;
k = 40;

X = randn(n,k)/sqrt(n);


s=svd(X,'econ');

a = 1-sqrt(k/n);
b = 1+sqrt(k/n);


figure(1),clf; hold on
[n1,x1]=hist(s, 10);
bar(x1,n1,'FaceColor','w','EdgeColor','b');
plot(s,zeros(k,1),'xb', 'MarkerSize',10);
plot(a,0,'xr', 'MarkerSize',20);
plot(b,0,'xr', 'MarkerSize',20);
axis([a-0.1, b+0.1, 0, max(n1)+1])
grid on; title(sprintf('n=%d, k=%d', n, k))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fast randomized svd

clear all; close all; rng(2017);

%% construnct true matrix

m = 1000;
n = 100;

k = 5;

eps1 = 0.01;
sA = [1:-0.1:1-0.1*(k-1), rand(1,n-k)*eps1]';
sA = sort(sA,'descend');

figure(11), clf;
plot(sA,'x-'); grid on;
title('singular value of A')

%
[u1,~,~]=svd(randn(m,n),'econ'); %random n orthonormal basis in R^m
[v1,~,~]=svd(randn(n,n),'econ'); %random n orthonormal basis in R^n

A = u1*diag(sA)*v1';
 
%% randomized svd: step 1, construnct Q
l = k*2+1;

Omega= randn(n,l);
Y = A* Omega;
[q1,r1,idx]=qr(Y,'vector');
Q= q1(:,1:l);

%% step 1, svd the reduced system
B = Q'*A;

[tilu,s,v]= svd(B,'econ');

u2 = Q*tilu;
s2 = diag(s);
v2 = v;

%% accuracy

A2 = u2*diag(s2)*v2';

diff_norm_op = norm(A-A2,2)
diff_bound = (1+4*sqrt((2*n)/(k-1)))*sA(k+1)

figure(12), clf; hold on;
plot(sA(1:l),'x-'); 
plot(s2,'o')
grid on;title('sA(1:l)')








