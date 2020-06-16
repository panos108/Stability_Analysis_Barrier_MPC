%%% Jury-Lee Theorem 2.3 paper TAC IEEE%%%
clc
clear
pp1=0;
mu=0.7;
[mini] = findminmax ([])

%%%%%%%%%%% for b = 0.5 ----> all r %%%%%%%%
%%%%%%%%%%% for r = 0.001 ----> b = 0.701 %%%%%%%
for par1 =.001
    k=par1;
    pp1=pp1+1;
    pp2=0;
    for par2 = 0.709%:0.001:1
        pp2=pp2+1;
        %%% Define system %%%
        z=tf('z');
        %%% Gain for plant %%%
        %   k=4.8766; %%%example 3 paper
        % k=6.7953;  %%% example diagonal S
        % % k=12.996;
        %    G=0.1*z/(z^2-1.8*z+0.81)*[5 2;3 4]*1;
        %
        % %   G = [0.2/(z-0.98),-.2/(z-0.92);0.3/(z-0.97),0.1/(z-0.91)]
        % % k=0.15;%%k for MPC
        % %%% Plant %%%
        % %     G=(.1*z-0.05)/(z^3-1.6*z^2+.55*z+0.072)*[5 2 1;1 4 5;2 1 4]*1;
        %     [A,B,C,D] = ssdata(balreal(G,0.001));
        %
        A =[0.7   0.3
            0.8   0.01];
        
        B=[1;0];
        
        
        C=[1 1.5];
        D = 0;
        
        
        
        G=tf(ss(A,B,C,D,[]));
        %%% Define sizes %%%
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
        
        Ju=minreal((z*eye(nx)-(A-A*L*C))\eye(nx)*B);
        Jy=minreal((z*eye(nx)-(A-A*L*C))\eye(nx)*L);
        
        %%% MPC
        
        
        N=2;
        M=2;
        E=zeros(nu,N*nu);
        E(1:nu,1:nu)=eye(nu);
        
        n2u=nu;
        n2p=np;
        
        % r = diag([k*ones(1,nu)]); %Weights on input moves (i.e. ||u(t)-u(t-1)||)
        q = eye(np); %Weights on output deviation from setpoint
        Q = sparse(kron(q,eye(N)));
        R = diag([k*ones(1,M*nu)]);
        %         LL=dare(A,B,eye(size(A)),k*eye(size(B,2)));
        %         Q=daug(Q(1:(N-1)*np,1:(N-1)*np),C*LL*C');
        %
        [Lambda Phi]=largematrices(N,M,nu,np,nx,A,B,C,D);%%X=Lambda*x+Phi*U
        F= full([Phi'*Q*Lambda]);
        %         Jy=1;
        %         Ju=0;
        %
        Gold=G;
        G=[daug(eye(np)*par2,-F)]*[zeros(np,np),G;Jy*eye(np),Ju+Jy*G]*[daug(eye(np)*par2,E)];
        [sys,g] = balreal(G);
        elim = (g<1e-8)         % Small entries of g are negligible states
        rsys = modred(sys,elim) % Remove negligible states
        
                G = rsys;

        [A,B,C,D]=ssdata(G);
        %%% Define sizes %%%
        nx=size(A,1);
        nu=size(B,2);
        np=size(C,1);
        nd=size(D,2);
        
        I=eye(nx);
        
        S= (full(Phi'*Q*Phi) + R+mu*mini*eye(size(R)))^-1;%%Hessian
        
        
        
        % S=k;
        K=S;
        
        
        
        
        
        ns=size(S,1);
        
        Phi=[eye(np),zeros(np);[eye(np)]*(1-1/z),zeros(np);zeros(np),eye(np);zeros(np),[eye(np)]*(1-1/z)];
        
        Ga=Phi*[G;eye(np)];
        [sysa,g] = balreal(Ga);
        elim = (g<1e-12) ;        % Small entries of g are negligible states
        rsys = modred(sysa,elim) ;% Remove negligible states
        [Aa,Ba,Ca,Da]=ssdata(rsys);%,0.0000001));
        na=size(Aa,1);
        nc=size(Ca,1);
        
        K2=diag(sdpvar(ns,1));
        K1=diag(sdpvar(ns,1));
        Ko=diag(sdpvar(ns,1));
        % %
        
        %         k1=1e-8;%sdpvar(1);
        %         K1=k1*eye(ns,ns);
        %
        %         k2=1e-8;%sdpvar(1);
        %         K2=	k2*eye(ns,ns);
        
        %         ko=sdpvar(1);
        %         Ko=ko*eye(ns,ns);
        %
        r1=sdpvar(1);
        R1=r1*eye(n2p,n2p);
        
        
        ro=r1;%sdpvar(1);
        Ro=ro*eye(n2p,n2p);
        
        X=[eye(na),zeros(na,np);Aa,Ba;Ca,Da];
        M11=[Ro,zeros(n2p,ns);zeros(ns,n2p),zeros(ns)];
        
        M13=[zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),Ko];
        
        M14=[zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),K2];
        
        M23=[zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),K1];
        
        M31=M13';
        
        M32=M23';
        
        M41=M14';
        
        M33=[-R1,zeros(n2p,ns);zeros(ns,n2p),-Ko*S^-1+(-Ko*S^-1)'];
        
        M34=[zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),-(K1+K2)*S^-1];
        
        M43=M34';
        
        MM=[M11,zeros(np),M13,M14;...
            zeros(np),zeros(np),M23,zeros(np);...
            M31,M32,M33,M34;...
            M41,zeros(np),M43,zeros(np)];
        
        % M2=[Ro,zeros(n2p),zeros(n2p),zeros(n2p);zeros(n2p),zeros(n2p),zeros(n2p),zeros(n2p);...
        %     1/0.1^2*eye(n2p),zeros(n2p),zeros(n2p),zeros(n2p);zeros(n2p),zeros(n2p),zeros(n2p),zeros(n2p)];
        %
        % kron(M2,[zeros(4,40),zeros(4);zeros(4,40),eye(4)])+kron(M1,[eye(40),zeros(40);zeros(40),zeros(40)])
        FF=[Ko>=(1e-8)*eye(size(Ko)),...
            K1>=(1e-8)*eye(size(K1)),...
            K2>=(1e-8)*eye(size(K2)),...
            R1>=(1e-8)*eye(size(R1)),...
            Ro>=(1e-8)*eye(size(Ro))];
        P=sdpvar(na,na);
        %         FF=[];
        GGG=[-P,zeros(na),zeros(na,nc);zeros(na),P,zeros(na,nc);zeros(nc,na),zeros(nc,na),MM];
        
        t=sdpvar(1);
        FF=[FF,X'*GGG*X<=-t*eye(size(X,2)),t>=1e-8,P>=1e-8];
        opt = sdpsettings('solver', 'mosek');%,'verbose',0);
        [p1]=optimize(FF,t,opt);
        
        % H1=value(R3)+value(R1)*(1-1/z)+value(R2)*(1-z)
        %
        % H2=value(Ko)+value(K1)*(1-1/z)+value(K2)*(1-z)
        
        
        max(eig(value(X'*GGG*X)))
        % MM=[zeros(np),zeros(np),Ko,K2;zeros(np),zeros(np),K1,zeros(np);Ko,K1,-2*Ko*S^-1,(K1+K2)*S^-1;K2,zeros(np),S^-1*(K1+K2),zeros(np)];
        % Phi'*MM*Phi
        % KK=[G;eye(3)]'*Phi'*MM*Phi*[G;eye(3)];
        
        
        % R3
        P
        p1
        results(pp1,pp2)=max(eig(value(X'*GGG*X)));
        
        
        if max(eig(value(X'*GGG*X)))>0
            results1(pp1,1)=4;
            break
        else
            results1(pp1,1)=0;
        end
    end
    if max(eig(value(X'*GGG*X)))>0
        results1(pp1,1)=4;
        break
    else
        results1(pp1,1)=0;
    end
end
% P1
