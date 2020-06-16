%%% Jury-Lee Theorem 2.3 paper TAC IEEE%%%
clc
clear
pp1=0;
mu=0.8;
[mini] = findminmax ([])

%%%%%%%%%%% for b = 0.5^2 ----> all r %%%%%%%%
%%%%%%%%%%% for r = 0.1 ----> b = 0.4140^2 %%%%%%%
%%%%%%%%%%% for n = 10    ----> b = 0.501^2 %%%%the most I can get
for par1 =0.1%.098:-0.001:0.0001
    k=par1;
    pp1=pp1+1;
    pp2=0;
    for par2 = 1.13:0.001:10%sqrt(0.5112)%0.5%:0.001:1
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
        G=[daug(eye(np)*0,-par2*F)]*[zeros(np,np),G;Jy*eye(np),Ju+Jy*G]*[daug(eye(np)*0,E)];
        [sys,g] = balreal(G);
        elim = (g<1e-8);         % Small entries of g are negligible states
        rsys = modred(sys,elim); % Remove negligible states
        G = rsys;
        
        [A,B,C,D]=ssdata(G);
        %%% Define sizes %%%
        nx=size(A,1);
        nu=size(B,2);
        np=size(C,1);
        nd=size(D,2);
        
        I=eye(nx);
        
        S= (full(Phi'*Q*Phi) + R)^-1;%%Hessian
        
        
        
        % S=k;
        K=S;
        
        
        
        
        
        ns=size(S,1);
        
        n = 10;   %%%%%%Size of Multiplier
        %         Phi1 = zeros(2*(n+1),2);
        for i =1:n
            
            Phi1(np*i+1:np*i+np,1:np) = [eye(np)]*(1-1/z^i);
            Phi1(np*i+np*n+1+np:np*i+np*n+2*np,np+1:2*np) = [eye(np)]*(1-1/z^i);
        end
        Phi1(1:np,1:np) = [eye(np)];
        Phi1(n*np+1+np:n*np+2*np,np+1:2*np) = [eye(np)];
        
        Phi = Phi1;
        %         Phi = [eye(np),         zeros(np);
        %             [eye(np)]*(1-1/z),zeros(np);
        %             [eye(np)]*(1-1/z^2),zeros(np);
        %             zeros(np)        ,eye(np);
        %             zeros(np)        ,[eye(np)]*(1-1/z);
        %             zeros(np)        ,[eye(np)]*(1-1/z^2)];
        
        Ga=Phi*[G;eye(np)];
        [sysa,g] = balreal(Ga);
        elim = (g<1e-12);         % Small entries of g are negligible states
        rsys = modred(sysa,elim); % Remove negligible states
        [Aa,Ba,Ca,Da]=ssdata(rsys);%,0.0000001));
        na=size(Aa,1);
        nc=size(Ca,1);
        
        
        for i = 1:2*n
            KK{i} = diag(sdpvar(ns,1));
        end
        
        
%         K4 = diag(sdpvar(ns,1));
%         K3 = diag(sdpvar(ns,1));
%         K2 = diag(sdpvar(ns,1));
%         K1 = diag(sdpvar(ns,1));
         Ko = diag(sdpvar(ns,1));
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
        
        
        
        

        
        
        row = size(Ro,1)+ns;
        col = size(Ro,2)+ns;
        
        
        
        
        M11 = [Ro,zeros(n2p,ns);zeros(ns,n2p),zeros(ns)];
        
        
        %         for i = 1:n+1
        %             M341{1:row,col*(i-1)+1:i*col} = [zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),KK{i}];
        %         end
        
        
        
        
        
        M14 = [zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),Ko];
%         
%         M15 = [zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),K2];
%         
%         M24 = [zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),K1];
%         
%         M16 = [zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),K4];
%         
%         M34 = [zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),K3];
%         
%         M41 = M14';
%         
%         M42 = M24';
%         
%         M61 = M16';
%         
%         M43 = M34';
%         
%         M51 = M15';
        
        M44 = [-R1,zeros(n2p,ns);zeros(ns,n2p),-Ko*S^-1+(-Ko*S^-1)'];
        
%         M45=[zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),-(K1+K2)*S^-1];
%         
%         M46=[zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),-(K3+K4)*S^-1];
%         
        
        %%%%% UPPER RIGHT
        Murtop = [];
        

        for i = 2:2:n*2
            Murtop = [Murtop, [zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),KK{i}]];

        end
        
        
        Murleft = [];
        

        for i = 1:2:n*2
            Murleft =  [Murleft;[zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),KK{i}]];

        end
        
        Mur0 =  zeros(n*row,n*col);
        
        Mur = [M14,Murtop;Murleft,Mur0];
        
        
       
        %%%%%
        %%%% LOWER LEFT
        
        Mll = Mur';
        %%%%
        
        %%%%% UPPER LEFT

        
        
                Multop = [];
        

        for i = 2:2:n*2
            Multop = [Multop, [zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),zeros(size(Ko))]];

        end
        
        
        Mulleft = [];
        

        for i = 1:2:n*2
            Mulleft =  [Mulleft;[zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),zeros(size(Ko))]];

        end
        
        
        Mul0 =  zeros(n*row,n*col);


        
                Mul = [M11,Multop;Mulleft,Mul0];

        %%%%%
        
        %%%%% LOWER RIGHT
      
        
        Mlrtop = [];
        

        for i = 1:2:n*2-1
            Mlrtop = [Mlrtop, [zeros(n2p,n2p),zeros(n2p,ns);zeros(ns,n2p),-(KK{i}+KK{i+1})*S^-1]];
        end
        
        
        Mlrleft = Mlrtop';
        

        
        
        Mlr0 =  zeros(n*row,n*col);


        
         Mlr = [M44,Mlrtop;Mlrleft,Mlr0];

%%%%%
       MM = [Mul, Mur;
             Mll, Mlr];
        
%                 M54=M45';
%         
%                 M64=M46';
%         
%                 MM1=[M11      , zeros(np), zeros(np), M14 , M15      , M16      ;...
%                     zeros(np), zeros(np), zeros(np), M24 , zeros(np), zeros(np);...
%                     zeros(np), zeros(np), zeros(np), M34 , zeros(np), zeros(np);...
%                     M41      , M42      , M43      , M44 , M45      , M46      ;...
%                     M51      , zeros(np), zeros(np), M54 , zeros(np), zeros(np);...
%                     M61      , zeros(np), zeros(np), M64 , zeros(np), zeros(np)];
%         
%         
        

                

       
        % M2=[Ro,zeros(n2p),zeros(n2p),zeros(n2p);zeros(n2p),zeros(n2p),zeros(n2p),zeros(n2p);...
        %     1/0.1^2*eye(n2p),zeros(n2p),zeros(n2p),zeros(n2p);zeros(n2p),zeros(n2p),zeros(n2p),zeros(n2p)];
        %
        % kron(M2,[zeros(4,40),zeros(4);zeros(4,40),eye(4)])+kron(M1,[eye(40),zeros(40);zeros(40),zeros(40)])
        FF=[Ko>=(1e-8)*eye(size(Ko)),...
            R1>=(1e-8)*eye(size(R1)),...
            Ro>=(1e-8)*eye(size(Ro))];
        
        for i =1:2*n
        FF = [FF,KK{i}>=(1e-8)*eye(size(KK{i}))];
        end
        
        
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
