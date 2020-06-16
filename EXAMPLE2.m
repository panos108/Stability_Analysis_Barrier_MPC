clear all
close all
clc
global A B N W y E Phi Lambda H F pF G S d_relax mu
% mu=0.896
% load notmon
% load dd
% load LIMITS_STATES_MU_0_05
% close all;
% d=x(7:12);

% systf=(1-s)/(s^2+s+1);
% [A,B,C,D]=tf2ss([0 -1 1],[1 1 1]);
A =[0.7   0.3
    0.8   0.01];

B=[1;0];


C=[1 1.5];
D = 0;

z = tf('z');

%-----------------------------------------------------
% Construct model
%------------------------------------------------
%Specify system here in state space form with A B C D
Ts = 1.0; %Sample interval
% A=[0.9 0;0.9 0.9]; B=[0.1;0.2]; C=[1 1]; D=0;
%Specify system dimensions
n_in = size(B,2); %# inputs
n_out = size(C,1); %# outputs
n_states = size(A,1); %# states
%------------------------------------------------------

% Specify horizons
%------------------------------------------------------
M = 2; %Control Horizon
N = M; %Prediction Horizon
%------------------------------------------------------
% Initial Matrix Filling
%------------------------------------------------------

%------------------------------------------------
% Quadratic cost weighting Q (error) and R (control) generation.
%------------------------------------------------


%-----------------------------------------------
% Create an observer
%-----------------------------------------------
nx=size(A,1);
nu=size(B,2);
np=size(C,1);
nd=size(D,2);

%%%observer

Robs=eye(np)*1;
Qobs=eye(nx);
[Lt,~,~] = dlqr(A',C',Qobs,Robs);
%  [L,prec,message] = place(A',C',[0.1 0.2 0.3 0.5 0.15 0.25 0.8 0.4 0.6]);
L=Lt';

Ju = minreal((z*eye(nx)-(A-A*L*C))\eye(nx)*B);
Jy = minreal((z*eye(nx)-(A-A*L*C))\eye(nx)*L);



% MPC



E=zeros(nu,N*nu);
E(1:nu,1:nu)=eye(nu);

n2u=nu;
n2p=np;
k = 0.1;
% r = diag([k*ones(1,nu)]); %Weights on input moves (i.e. ||u(t)-u(t-1)||)
q = eye(np); %Weights on output deviation from setpoint
Q = sparse(kron(q,eye(N)));
R = diag([k*ones(1,M*nu)]);
%         LL=dare(A,B,eye(size(A)),k*eye(size(B,2)));
%         Q=daug(Q(1:(N-1)*np,1:(N-1)*np),C*LL*C');
%
[Lambda Phi]=largematrices(N,M,nu,np,nx,A,B,C,D);%%X=Lambda*x+Phi*U
F= full([Phi'*Q*Lambda]);H = full(Phi'*Q*Phi) + R;%%Hessian

%%%obj=0.5U'HU-Fx+mu*B(x,U)%%
W=[ones(N*n_in,1);.5*ones(N*n_in,1)];
% E(N*n_states+1,:)=rand(1,N*n_states); %%%%
opt=optimset('display','iter','MaxFunEvals',100000,'MaxIter',10000000);

I=eye(size(H,1));

Phi=full(Phi);

mu = 1.001;
d_relax(1:size(E*Phi)) = 0.00000001;
uqp = zeros(N*n_in,1);


G1=tf(ss(A,B,C,D,[]));




Gold=G1;
par2 = 0;0.25;
G1=[daug(eye(np)*par2,-F)]*[zeros(np,np),G1;Jy*eye(np),Ju+Jy*G1]*[daug(eye(np)*par2,[1 0])];
[sys,g] = balreal(G1);
elim = (g<1e-8);       % Small entries of g are negligible states
rsys = modred(sys,elim); % Remove negligible states


[A1,B1,C1,D1]=ssdata(rsys);

x = [1;-1];
xr = x+rand;
w = 0;
kk = 2.9;
y1 = C*xr+1;
rrr = rand;
for i = 1:1
    
x = [.1;-.1];
xr = x;randn;
mu = 0.4;mu - 0.2
m(i) = mu;
w = 1;
y1 = C*xr+1;
kk = kk;% + .1;
LU=eye(N*n_in);

    for j =1:200
        
        %             x(1,1)=x1(i);
        %             %     for j=1:length(xx)
        %             x(2,1)=x2(k);%x2(j);;
        y=-kk*F*x;
          [U,j1,h]=fmincon(@Barrier2,zeros(2,1),[],[],[],[],[],[],[],[]);
%        [U,FVAL,EXITFLAG,OUTPUT,lambda] = quadprog(H,-y',[LU;-LU],W,[],[],[],[],[],opt);
        u = U(1);
        
        xxr1(j,i) = xr(1);
        xxx1(j,i) = x(1);
        
        xxr2(j,i) = xr(2);
        xxx2(j,i) = x(2);
        
        
        
        y1  = C*xr+par2*w;
        xr = A*xr+B*u;
        x = (A-A*L*C)*x+B*u+A*L*(y1);
        %     xx = A1*xx+B1*[w;U];
        yy(j,i) =y1;
        %     x = yy(2:3,j);
        
        if norm(par2*y1)<=1
            w = par2*y1;
        elseif (par2*y1)>1
            w = (par2*y1-1);
        else
            w = (par2*y1+1);
        end
%         w = 0.25*rand;
        s(j) = norm(w)/norm(par2*y1);
        uu(j) = u;
    end
end
% plot(xxr');hold on;plot(xxx','--')
plot(yy)

% count = load('count.dat');
t = (1:size(yy,1))';
ym = mean(yy,2);
eym = std(yy,1,2);


ylo = ym - eym;
yhi = ym + eym;

t = (1:size(yy,1))';
xm1 = mean(xxr1,2);
exm1 = std(xxr1,1,2);


xlo1 = xm1 - exm1;
xhi1 = xm1 + exm1;

t = (1:size(yy,1))';
xm2 = mean(xxr2,2);
exm2 = std(xxr2,1,2);


xlo2 = xm2 - exm2;
xhi2 = xm2 + exm2;


hp = patch([t; t(end:-1:1); t(1)], [ylo; yhi(end:-1:1); ylo(1)], 'r');
hold on;
h1 = line(t,ym);

set(hp, 'facecolor', [1 0.8 0.8], 'edgecolor', 'none','facealpha',[0.3]);
set(h1, 'color', 'r', 'marker', 'x','linewidth',3);


% hp = patch([t; t(end:-1:1); t(1)], [xlo1; xhi1(end:-1:1); xlo1(1)], 'r');
% hold on;
% h3 = line(t,xm1);
% 
% set(hp, 'facecolor', [.8 0.8 1], 'edgecolor', 'none','facealpha',[0.3]);
% set(h3, 'color', 'b', 'marker', '*','linewidth',3);
% 
% hp = patch([t; t(end:-1:1); t(1)], [xlo2; xhi2(end:-1:1); xhi2(1)], 'r');
% hold on;
% h2 = line(t,xm2);
% 
% set(hp, 'facecolor', [0.8 1 0.8], 'edgecolor', 'none','facealpha',[0.3]);
% set(h2, 'color', 'g', 'marker', '.','linewidth',3);
% 


% 

legend([h1 h3 h2],{'\it{Output}','\it{x_1}','\it{x_2}'});  % Only the blue and green lines appear


