function [mini] = findminmax (mo)
% load dd
global b
global A B N G 

%%%%------------------------------------------
%%% Create the matrices of the system
% A=[.5,.1;-.1,1];
% B=[0.05 ,-0.02;.1,0.01];
%  C=[1 0;0 1];
% D=[0];
% load lol

k=mo;
% A =[0.16   0.31
%     0.79   0.53];
A =[0.7   0.3
    0.8   0.01];

B=[1;0];


C=[1 1.5];
D = 0;




%-----------------------------------------------------
% Construct model
%------------------------------------------------
%Specify system here in state space form with A B C D
Ts = 1.0; %Sample interval
%Specify system dimensions
n_in = size(B,2); %# inputs
n_out = size(C,1); %# outputs
n_states = size(A,1); %# states
%------------------------------------------------------

% Specify horizons
%------------------------------------------------------
M =2; %Control Horizon
N = M; %Prediction Horizon
%------------------------------------------------------


%%%%----create the matrices and the constraints of the problem--------
%-----------------------------------------------

opt=optimoptions(@fmincon,'Algorithm','interior-point','Display','off','MaxFunctionEvaluations',1000000000,'MaxIterations',1500000000);



opts = optimoptions('fmincon','Display','off','Algorithm','interior-point');

global LU
LU=eye(N*n_in);

x=rand(n_states*N+n_in*N,1000);


b = [ones(size(LU,2),1)-0;.5*ones(size(LU,2),1)];
G = [LU;-LU];
[xv2,hl,j,~,l]=fmincon(@objectived_min,rand(1,size(G,2)),...
    [G],b,[],[],[],[],[],opts);%
xx(:,1)=xv2;

OPTT.InitialPopulationMatrix=xx;
%
OPTT.PlotFcn = {@gaplotbestf,@gaplotbestindiv};

problem = createOptimProblem('fmincon','x0',rand(1,size(G,2)),...
    'objective',@objectived_min,'Aineq',[G],'bineq',[b],'options',opts);

gs = MultiStart;
[xmin,f] = run(gs,problem,200);


mini = f/2;
end