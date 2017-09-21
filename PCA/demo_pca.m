% loss of consistency for large-p PCA

clear all; close all; rng(2017);

%% the signal + noise case

n=2000;
%n=4000;
%n=1000;
%p=100;
p=1000;

u=zeros(p,1);
u(1)=1;


sigz=.5;

SNR=1/sigz^2

SNR_th=sqrt(p/n)

%%

x=u*randn(1,n);
y=x + sigz*randn(p,n);

Sy = (y*y')/n;

%
[v,d]=eig(Sy);
[lambda_y, tmp]=sort(diag(d),'descend');
psi_y=v(:,tmp);

lambda_y_pop = norm(u)^2+sigz^2

figure(2),clf;
%hist(lambda_y,100);
[nout,xout]=hist(lambda_y,100); 
bar(xout,nout,'Edgecolor','b','Facecolor','w'); grid on;
hold on; 
plot(lambda_y_pop,0,'xr');
plot(sigz^2,0,'xr');
grid on; title(sprintf('lambda y, n=%d, p=%d', n, p));
axis([0,2,0,50])

%
figure(3),clf;hold on;
p1=psi_y(:,1);
if p1'*u < 0
    p1=-p1;
end
plot(p1,'xb')
plot(u,'.-r')
title('first eigenvetor of Sy')
grid on;

%%

%%


%% the null case

n=2000;
n=5000;

%n=10000;
p=1000;

x=randn(p,n);
Sx= (x*x')/n;

%
[v,d]=eig(Sx);
[lambda_x, tmp]=sort(diag(d),'descend');
psi_x=v(:,tmp);

lambda_x_pop = 1

%%

gm=p/n;
a=(1-sqrt(gm))^2;
b=(1+sqrt(gm))^2;
tt=max(0,a-0.1):1e-3:b+0.1;

func_MP=@(t) ( (t<b).*(t>a) ).*sqrt((b-t).*(t-a))./(2*pi*gm*t);

pp=func_MP(tt);
%%

figure(4),clf;
%hist(lambda_y,100);
[nout,xout]=hist(lambda_x,100); 
bar(xout,nout/p,'Edgecolor','b','Facecolor','w'); grid on;
hold on; 
plot(lambda_x_pop,0,'xr');
plot(tt,pp*(xout(2)-xout(1)),'--b','Linewidth',2);
grid on; title(sprintf('lambda x, null case n=%d, p=%d', n, p));
