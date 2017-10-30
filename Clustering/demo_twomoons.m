% this demo shows the two circle example, compare k-means and spectral
% clustering

clear all; close all; rng(2017);

%% data 

n=200;

dim=2;

gapx=0.05; %0.1;

epsx=0.05; %0.1;

hh=1/(n/2);
tt=(hh/2:hh:1)';
x1=[cos(tt*pi),sin(tt*pi)];
x1(:,1)=x1(:,1)-.5;
x1(:,2)=x1(:,2)-gapx;

x2=[cos(tt*pi),sin(tt*pi)];
x2=-x2;
x2(:,1)=x2(:,1)+.5;
x2(:,2)=x2(:,2)+gapx;

x=cat(1,x1,x2);

x=x+randn(size(x))*epsx;

y=[ones(n/2,1);ones(n/2,1)*2];

%
figure(1),clf; 
scatter(x(:,1),x(:,2),80,y); colormap(jet);
grid on; title('data (color by labels)')

%%

%% spectral clustering

dis=pdist(x);

sig=0.25; 
W= exp(-squareform(dis.^2)/(2*sig^2));

figure(2),clf;
imagesc(W);colorbar();
title('W')

%% L_rw

dW=sum(W,2);

[v,d]=eig(W,diag(dW));
[lambda1,tmp]=sort(diag(d),'descend');
psi1=v(:,tmp);

%

figure(3),clf;
scatter(x(:,1),x(:,2),80,psi1(:,2));
grid on;colorbar();
title('psi2')


figure(5),clf;
bar(lambda1(1:10));grid on;
title('lambda')


%%

figure(4),clf;
scatter(psi1(:,2),psi1(:,3),80,y);
grid on; colormap(jet);
title('spectral embedding (color by labels)')



%% look at the whole spectrum

figure(6),clf;
imagesc(psi1(:,1:10));
title('psi')

figure(7),clf;
scatter3(psi1(:,2),psi1(:,3),psi1(:,4),80,y);
grid on; colormap(jet);
title('spectral embedding (color by labels)')

%%
%% k-means

numcluster=2;

opts_kmeans=statset('Display','final');%'iter');

%y2 = kmeans(x,numcluster,'Replicates',1,'Options',opts_kmeans);
y2 = kmeans(x,numcluster,'Replicates',10,'Options',opts_kmeans);


figure(8),clf;
scatter(x(:,1),x(:,2),80,y2);
grid on;colormap(jet);
title('cluster by kmeans')