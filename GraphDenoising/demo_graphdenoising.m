
clear all; close all; rng(2017)

%% smooth signal on a graph

n=100;
hh=1/n;

tt = (hh/2:hh:1)';

% signal
f = exp(-(tt - 0.4).^2/(2*0.1^2)) + .5* exp(-(tt - 0.75).^2/(2*0.08^2)) ;

% add noise
sigx = 0.1;
x = f + randn(n,1)*sigx;

figure(1),clf; hold on;
plot(tt, f,'.-');
plot(tt, x,'o');
grid on; title('f and x')

%% estimation

% just use observation
f1=x;

% graph denoising
dis = pdist(tt);
W = exp(- squareform(dis).^2/(2*0.05^2));

figure(2),clf;
imagesc(W); title('W')

dW=sum(W,2);
P = diag(dW)\W;

f2 = P*x;

figure(3),clf; hold on;
plot(tt, f,'.-');
plot(tt, x,'o');
plot(tt, f2,'x-');
grid on; title('f, x, fest')

% mse
mse1 = sum( (f1-f).^2)/n
mse2 = sum( (f2-f).^2)/n

%% 


%% manifold plus noise data

% non-local means applied to toy-model

n=100;
hh=1/n;

tt = (hh/2:hh:1)';

% 
f = [cos(tt*pi), sin(tt*pi)];

sigx = 0.05;
x = f + randn(size(f))*sigx;



figure(11),clf; hold on;
scatter(f(:,1),f(:,2),'.b');
scatter(x(:,1),x(:,2),'og');
grid on; title('f and x')

%% graph denoising

sig = 0.1; 

dis = pdist(x);
W = exp(- squareform(dis).^2/(2*sig^2));

figure(12),clf;
imagesc(W); title('W')

%%
dW=sum(W,2);
P = diag(dW)\W;

f2 = P*x;

figure(13),clf; hold on;
scatter(f(:,1),f(:,2),'.b');
scatter(x(:,1),x(:,2),'og');
scatter(f2(:,1),f2(:,2),'xr');
grid on; title('f, x, fest')
