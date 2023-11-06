%% Stochastic square lattice forest-grass-fire-ash model 
%% parameters, adjacency and initialisation   
L=100; N=L^2; 
pg=.9;pf=.1; %fire spreading prob in grass and forest
    % p = ρ/(μ + ρ) => ρ = μ  p/(1-p)

lambda=5.;        %ash to grass conversion rate
mu=10^6;          %fire to ash conversion rate 
rhog=mu*pg/(1-pg);      %fire contagion rate in grass
rhof=mu*pf/(1-pf);      %fire contagion rate in forest
alpha=0.03;      %facilitated forest growth rate
beta=0.0002;      %spontaneous forest growth rate
phi=.1/N;      %spontaneous ignition rate
gamma=0.02;       %spontaneous forest mortality 
fG0=1.;nG0=N*fG0; %initial grass fraction of domain (rest is forest)

%construct rate matrix/tensor 
    %M1: spontaneous conversion rates (i->j)
    %M2: facilitated conversion rates (ik->jk)
%   1 2 3 4
%X=[F,G,B,A]
M1=[0 gamma 0 0; ...
    beta 0 phi 0;...
    0 0 0 mu; ...
    beta lambda 0 0];
M2=zeros(4,4,4);
M2(2,3,3)=rhog;  %GB->BB
M2(1,3,3)=rhof;  %FB->BB
M2(4,1,1)=alpha; %FA->FF
M2(2,1,1)=alpha; %FG->FF

rates.M1=M1;
rates.M2=M2;

%max simulation time/iteration
tmax=500;
imax=2*10^7;

%adjacency matrix for square lattice
ii = speye(L,L);
dd = spdiags([ones(L,1),ones(L,1)],[-1 1],L,L);
Adj = logical(kron(dd,ii)+kron(ii,dd)); %use boolean for performance
%add links for periodic BC
left=1:L; right=(1:L)+L*(L-1); top=(L:L:L^2)-L+1;bottom=(L:L:L^2);
bidx = sub2ind(size(Adj), [left,top,right,bottom],[right,bottom,left,top]);
Adj(bidx)=true;

%initialise F,G,B,A
B=false(N,1);
G=[true(nG0,1);false(N-nG0,1)];
G=G(randperm(N));
A=false(N,1);
F=~G;
X0=sparse([F,G,B,A]);

%plot colours
clrF=[34. 71. 42]./255;
clrG=[189. 208. 156.]./255;
clrB=[217. 145. 22.]./255;
clrA=[178. 178. 178.]./255;
clr={clrF;clrG;clrB;clrA};%colours;

%% regular sim

plt.di=1000;
plt.clr=clr;
plt.vid=[];

tmax=300;
dat=KMC(X0,Adj,rates,L,tmax,imax,plt,0);

%% reconstruct video from data

nvars=4;
nt=size(dat,1)-N;
figure;
X=dat(1:N,4);
t=0;t0=0;tend=300;
dT=1;
timeelapsed=0;
cmap=[0 0 0; cell2mat(plt.clr)];

for i=1:nt
    if timeelapsed>=dT
        imagesc(reshape(X,L,L));
        colormap(gca,cmap);caxis([0,nvars]);
        axis square; axis off;
        box on;
        title(sprintf('t=%.6f/%u',t,tmax));
        drawnow;%waitforbuttonpress;%
        t0=t;
    end
        X(dat(N+i,2))=dat(N+i,4);
        t=dat(N+i,1);
        timeelapsed=t-t0;
        if t>=tend
            break
        end
end