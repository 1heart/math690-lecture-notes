
clear all; close all; 
rng(201709);

%% swiss roll data
N=2000;

  tt = (3*pi/2)*(1+2*rand(1,N));  
  tt=sort(tt);
  height = 21*rand(1,N);
  X = [tt.*cos(tt); height; tt.*sin(tt)];

% 
figure(1),clf;
  scatter3(X(1,:),X(2,:),X(3,:),12,tt,'+');
  view([12 20]); grid on; 
 drawnow;
 
%% PCA
x=X';
mx=mean(x,1);
xx=bsxfun(@minus, x, mx);

[ux,sx,vx]=svd(xx/sqrt(N),'econ');

diag(sx)
 
y_PCA=ux(:,1:2);

%% MDS
D=squareform(pdist(X'));

figure(19),clf;
imagesc(D);colorbar();
title('Euclidean Distance')
%%

M = -.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2));
M = (M + M')/2;

ndims=10;
%[v,d]=eigs(M, ndims, 'LR');
[v,d]=eig(M);

[lambda,tmp]=sort(diag(d),'descend');
y_MDS=v(:,tmp(1:3));

 %%
figure(20), clf;
scatter(y_PCA(:,1),y_PCA(:,2),40,tt,'+')
grid on; title('PCA')
 
figure(21), clf;
scatter(y_MDS(:,1),y_MDS(:,2),40,tt,'+')
grid on; title('MDS')
 
%% Iso map

%% construnct nn graph
kNN=50;
%kNN=30;
%kNN=20;
%kNN=10; %what is the effect of setting knn?


[idx,dis]=knnsearch(X',X','k',kNN);

I=repmat((1:N)',[1,kNN]);
I=I(:);
J=idx(:);
%V=exp(-dis.^2/(2*sig2g));
V=ones(size(J));
W=sparse(I,J,V,N,N);

%A=max(W,W');
A=min(W,W');

%%  Compute shortest paths
G = graph(A,'upper');
D= distances(G);

figure(2),clf;
imagesc(A);title('A')

figure(3),clf;
imagesc(D);title('Graph Shortest Path Distance')

%% classical MDS to isomap

M = -.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2));
M = (M + M')/2;

ndims=10;
%[v,d]=eigs(M, ndims, 'LR');
[v,d]=eig(M);

[lambda,tmp]=sort(diag(d),'descend');
Y=v(:,tmp(1:3));

figure(10), clf;
scatter(Y(:,1),Y(:,2),40,tt,'+')
grid on; title('iso map') %look at fingers at the boundary

figure(11),clf;
plot(lambda(1:20),'x-');
grid on; title('lambda mds');

figure(12), clf;
scatter3(Y(:,1),Y(:,2),Y(:,3),40,tt,'+')
grid on; title('iso map 3')

return;

%% laplacian eigenmap 

%%
%W=full(A);

sig2g=2^2;

kNN=500;
[idx,dis]=knnsearch(X',X','k',kNN);

I=repmat((1:N)',[1,kNN]);
I=I(:);
J=idx(:);
V=exp(-dis.^2/(2*sig2g));

W=sparse(I,J,V,N,N);
W=(W+W')/2;
W=full(W);

figure(32),clf;
imagesc(W);colorbar();

figure(33),clf;
scatter3(X(1,:),X(2,:),X(3,:),40, W(:,500),'o','filled')



%%
dW=sum(W,2);


tic,
[v,d]=eigs(W,diag(dW),10,'la');
toc

[lambda,tmp]=sort(diag(d),'descend');
Y_Lap=v(:,tmp);


figure(31), clf;
scatter3(Y_Lap(:,2),Y_Lap(:,5),Y_Lap(:,6),40,tt,'+')
%scatter3(Y_Lap(:,5),Y_Lap(:,6),Y_Lap(:,7),40,tt,'+')
grid on; title('Laplace Eigen Map') %look at fingers at the boundary


