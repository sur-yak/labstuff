function [x,y] = SolutionLambdaOmegaCubic_NetworkLinear(t,params,Fx,Fy,etax,etay,NSETYPE,NUMMETH,initcond,Mod2,FLdaOmegaCub,GLdaOmegaCub)


% NSETYPE:
            %  0: No noise
            %  1: Piecewise-constant noise (Mean=0, Var=Dx, Var=Dy)
            %  2: White (Gaussian) noise (Honeycut, 1992; Only for RK2)
            %  3: Ornstein-Uhlenbeck noise


%NUMMETH = 2;
            % 1: Runge-Kutta Order 2
            % 2: Runge-Kutta Order 4

 % Time-step
 
dt = t(2)-t(1);

Nosc = params(1);
lda = params(2:Nosc+1);
b = params(Nosc+2:2*Nosc+1);
omega = params(2*Nosc+2:3*Nosc+1);
a = params(3*Nosc+2:4*Nosc+1);
ldah = params(4*Nosc+2:5*Nosc+1);
bh = params(5*Nosc+2:6*Nosc+1);
omegah = params(6*Nosc+2:7*Nosc+1);
ah = params(7*Nosc+2:8*Nosc+1);
Dx = params(8*Nosc+2:9*Nosc+1);
Dy = params(9*Nosc+2:10*Nosc+1);
Tdur = params(10*Nosc+2);
alphavec = params(10*Nosc+3:10*Nosc+3+Nosc*Nosc-1);
betavec = params(10*Nosc+3+Nosc*Nosc:10*Nosc+3+2*Nosc*Nosc-1);
alphaimagvec = params(10*Nosc+3+2*Nosc*Nosc:10*Nosc+3+3*Nosc*Nosc-1);
betaimagvec = params(10*Nosc+3+3*Nosc*Nosc:10*Nosc+3+4*Nosc*Nosc-1);
ton = params(10*Nosc+3+4*Nosc*Nosc:11*Nosc+3+4*Nosc*Nosc-1);
tconn = params(11*Nosc+3+4*Nosc*Nosc);
% ton = params(10*Nosc+3+Nosc*Nosc:11*Nosc+3+Nosc*Nosc-1);
% tconn = params(11*Nosc+3+Nosc*Nosc);
xovec = initcond(1:Nosc);
yovec = initcond(Nosc+1:2*Nosc);

alphaconn = reshape(alphavec,Nosc,Nosc);
alpha0 = zeros(Nosc);
alpha = alpha0;
betaconn = reshape(betavec,Nosc,Nosc);
beta0 = zeros(Nosc);
beta = beta0;
alphaimagconn = reshape(alphaimagvec,Nosc,Nosc);
alphaimag0 = zeros(Nosc);
alphaimag = alphaimag0;
betaimagconn = reshape(betaimagvec,Nosc,Nosc);
betaimag0 = zeros(Nosc);
betaimag = betaimag0;



Tdur;

% Piecewise-constant Input functions

if NSETYPE == 1
    for k=1:Nosc
        Fx(k,:) = sqrt(Dx(k))*Fx(k,:);
        Fy(k,:) = sqrt(Dy(k))*Fy(k,:);
    end
end

% Numerical simulations

x = zeros(Nosc,length(t));
y = zeros(Nosc,length(t));

x(:,1) = xovec;
y(:,1) = yovec;

if NUMMETH == 1
    kx1 = zeros(1,Nosc);
    ky1 = zeros(1,Nosc);
    kx2 = zeros(1,Nosc);
    ky2 = zeros(1,Nosc);
    ax = zeros(1,Nosc);
    ay = zeros(1,Nosc);
    if NSETYPE == 0
        for j=1:length(t)-1
            for k=1:Nosc
                if t(j) < ton(k)
                    x(k,j) = x(k,1);
                    y(k,j) = y(k,1);
                end
            end
            if j == floor(tconn/dt)
                alpha = alphaconn;
                beta = betaconn;
                alphaimag = alphaimagconn;
                betaimag = betaimagconn;
            end
            for k=1:Nosc
                mod = Mod2(x(k,j),y(k,j));
                kx1(k) = FLdaOmegaCub(x(k,j),y(k,j),mod,lda(k),b(k),omega(k),a(k));
                kx1(k) = kx1(k)+sum(alpha(k,:).*x(:,j)');
                kx1(k) = kx1(k)-sum(alphaimag(k,:).*y(:,j)');
                ky1(k) = GLdaOmegaCub(x(k,j),y(k,j),mod,ldah(k),bh(k),omegah(k),ah(k)); 
                ky1(k) = ky1(k)+sum(beta(k,:).*y(:,j)');
                ky1(k) = ky1(k)+sum(betaimag(k,:).*x(:,j)');
            end
            for k=1:Nosc
                ax(k) = x(k,j)+kx1(k)*dt;    
                ay(k) = y(k,j)+ky1(k)*dt;  
            end
            for k=1:Nosc
                mod = Mod2(ax(k),ay(k));
                kx2(k) = FLdaOmegaCub(ax(k),ay(k),mod,lda(k),b(k),omega(k),a(k));
                kx2(k) = kx2(k)+sum(alpha(k,:).*ax(:)');
                kx2(k) = kx2(k)-sum(alphaimag(k,:).*ay(:)');
                ky2(k) = GLdaOmegaCub(ax(k),ay(k),mod,ldah(k),bh(k),omegah(k),ah(k));
                ky2(k) = ky2(k)+sum(beta(k,:).*ay(:)');
                ky2(k) = ky2(k)+sum(betaimag(k,:).*ax(:)');
            end
            for k=1:Nosc
                x(k,j+1) = x(k,j)+(kx1(k)+kx2(k))*dt/2;
                y(k,j+1) = y(k,j)+(ky1(k)+ky2(k))*dt/2;
            end
        end
    elseif NSETYPE == 1  
        for j=1:length(t)-1
            for k=1:Nosc
                if t(j) < ton(k)
                    x(k,j) = x(k,1);
                    y(k,j) = y(k,1);
                end
            end
            if j == floor(tconn/dt)
                alpha = alphaconn;
                beta = betaconn;
                alphaimag = alphaimagconn;
                betaimag = betaimagconn;
            end
            for k=1:Nosc
                mod = Mod2(x(k,j),y(k,j));
                kx1(k) = FLdaOmegaCub(x(k,j),y(k,j),mod,lda(k),b(k),omega(k),a(k))+Fx(k,j);
                kx1(k) = kx1(k)+sum(alpha(k,:).*x(:,j)');
                kx1(k) = kx1(k)-sum(alphaimag(k,:).*y(:,j)');
                ky1(k) = GLdaOmegaCub(x(k,j),y(k,j),mod,ldah(k),bh(k),omegah(k),ah(k))+Fy(j); 
                ky1(k) = ky1(k)+sum(beta(k,:).*y(:,j)');
                ky1(k) = ky1(k)+sum(betaimag(k,:).*x(:,j)');
            end
            for k=1:Nosc
                ax(k) = x(k,j)+kx1(k)*dt;    
                ay(k) = y(k,j)+ky1(k)*dt;  
            end
            for k=1:Nosc
                mod = Mod2(ax(k),ay(k));
                kx2(k) = FLdaOmegaCub(ax(k),ay(k),mod,lda(k),b(k),omega(k),a(k))+Fx(k,j+1);
                kx2(k) = kx2(k)+sum(alpha(k,:).*ax(:)');
                kx2(k) = kx2(k)-sum(alphaimag(k,:).*ay(:)');
                ky2(k) = GLdaOmegaCub(ax(k),ay(k),mod,ldah(k),bh(k),omegah(k),ah(k))+Fy(j+1);
                ky2(k) = ky2(k)+sum(beta(k,:).*ay(:)');
                ky2(k) = ky2(k)+sum(betaimag(k,:).*ax(:)');
            end
            for k=1:Nosc
                x(k,j+1) = x(k,j)+(kx1(k)+kx2(k))*dt/2;
                y(k,j+1) = y(k,j)+(ky1(k)+ky2(k))*dt/2;                
            end
        end
    elseif NSETYPE == 2
        for j=1:length(t)-1
            for k=1:Nosc
                if t(j) < ton(k)
                    x(k,j) = x(k,1);
                    y(k,j) = y(k,1);
                end
            end
            if j == floor(tconn/dt)
                alpha = alphaconn;
                beta = betaconn;
                alphaimag = alphaimagconn;
                betaimag = betaimagconn;
            end
            for k=1:Nosc
                mod = Mod2(x(k,j),y(k,j));
                kx1(k) = FLdaOmegaCub(x(k,j),y(k,j),mod,lda(k),b(k),omega(k),a(k));
                kx1(k) = kx1(k)+sum(alpha(k,:).*x(:,j)');
                kx1(k) = kx1(k)-sum(alphaimag(k,:).*y(:,j)');
                ky1(k) = GLdaOmegaCub(x(k,j),y(k,j),mod,ldah(k),bh(k),omegah(k),ah(k)); 
                ky1(k) = ky1(k)+sum(beta(k,:).*y(:,j)');
                ky1(k) = ky1(k)+sum(betaimag(k,:).*x(:,j)');
            end
            for k=1:Nosc
                ax(k) = x(k,j)+kx1(k)*dt;    
                ay(k) = y(k,j)+ky1(k)*dt;  
                ax(k) = ax(k)+sqrt(2*Dx(k)*dt)*etax(k,j);  
                ay(k) = ay(k)+sqrt(2*Dy(k)*dt)*etay(k,j); 
            end
            for k=1:Nosc
                mod = Mod2(ax(k),ay(k));
                kx2(k) = FLdaOmegaCub(ax(k),ay(k),mod,lda(k),b(k),omega(k),a(k));
                kx2(k) = kx2(k)+sum(alpha(k,:).*ax(:)');
                kx2(k) = kx2(k)-sum(alphaimag(k,:).*ay(:)');
                ky2(k) = GLdaOmegaCub(ax(k),ay(k),mod,ldah(k),bh(k),omegah(k),ah(k));
                ky2(k) = ky2(k)+sum(beta(k,:).*ay(:)');
                ky2(k) = ky2(k)+sum(betaimag(k,:).*ax(:)');
            end
            for k=1:Nosc
                x(k,j+1) = x(k,j)+(kx1(k)+kx2(k))*dt/2;
                y(k,j+1) = y(k,j)+(ky1(k)+ky2(k))*dt/2;                
                x(k,j+1) = x(k,j+1)+sqrt(2*Dx(k)*dt)*etax(k,j);
                y(k,j+1) = y(k,j+1)+sqrt(2*Dy(k)*dt)*etay(k,j); 
            end
        end
    end    
elseif NUMMETH == 2
    kx1 = zeros(1,Nosc);
    ky1 = zeros(1,Nosc);
    kx2 = zeros(1,Nosc);
    ky2 = zeros(1,Nosc);
    kx3 = zeros(1,Nosc);
    ky3 = zeros(1,Nosc);
    kx4 = zeros(1,Nosc);
    ky4 = zeros(1,Nosc);
    ax1 = zeros(1,Nosc);
    ay1 = zeros(1,Nosc);
    ax2 = zeros(1,Nosc);
    ay2 = zeros(1,Nosc);
    ax3 = zeros(1,Nosc);
    ay3 = zeros(1,Nosc);
    if NSETYPE == 0
        for j=1:length(t)-1
            for k=1:Nosc
                if t(j) < ton(k)
                    x(k,j) = x(k,1);
                    y(k,j) = y(k,1);                 
                end
            end          
            if j == floor(tconn/dt)
                alpha = alphaconn;
                beta = betaconn;
                alphaimag = alphaimagconn;
                betaimag = betaimagconn;   
            end
            for k=1:Nosc
                mod = Mod2(x(k,j),y(k,j));
                kx1(k) = FLdaOmegaCub(x(k,j),y(k,j),mod,lda(k),b(k),omega(k),a(k));                
                kx1(k) = kx1(k)+sum(alpha(k,:).*x(:,j)');
                kx1(k) = kx1(k)-sum(alphaimag(k,:).*y(:,j)');
                ky1(k) = GLdaOmegaCub(x(k,j),y(k,j),mod,ldah(k),bh(k),omegah(k),ah(k));         
                ky1(k) = ky1(k)+sum(beta(k,:).*y(:,j)');
                ky1(k) = ky1(k)+sum(betaimag(k,:).*x(:,j)');
            end           
             for k=1:Nosc
                ax1(k) = x(k,j)+kx1(k)*dt/2;    
                ay1(k) = y(k,j)+ky1(k)*dt/2;  
            end
            for k=1:Nosc
                mod = Mod2(ax1(k),ay1(k));
                kx2(k) = FLdaOmegaCub(ax1(k),ay1(k),mod,lda(k),b(k),omega(k),a(k));
                kx2(k) = kx2(k)+sum(alpha(k,:).*ax1(:)');
                kx2(k) = kx2(k)-sum(alphaimag(k,:).*ay1(:)');
                ky2(k) = GLdaOmegaCub(ax1(k),ay1(k),mod,ldah(k),bh(k),omegah(k),ah(k));
                ky2(k) = ky2(k)+sum(beta(k,:).*ay1(:)');
                ky2(k) = ky2(k)+sum(betaimag(k,:).*ax1(:)');
            end
            for k=1:Nosc
                ax2(k) = x(k,j)+kx2(k)*dt/2;    
                ay2(k) = y(k,j)+ky2(k)*dt/2;  
            end
            for k=1:Nosc
                mod = Mod2(ax2(k),ay2(k));
                kx3(k) = FLdaOmegaCub(ax2(k),ay2(k),mod,lda(k),b(k),omega(k),a(k));
                kx3(k) = kx3(k)+sum(alpha(k,:).*ax2(:)');
                kx3(k) = kx3(k)-sum(alphaimag(k,:).*ay2(:)');
                ky3(k) = GLdaOmegaCub(ax2(k),ay2(k),mod,ldah(k),bh(k),omegah(k),ah(k));
                ky3(k) = ky3(k)+sum(beta(k,:).*ay2(:)');
                ky3(k) = ky3(k)+sum(betaimag(k,:).*ax2(:)');
            end
            for k=1:Nosc
                ax3(k) = x(k,j)+kx3(k)*dt;    
                ay3(k) = y(k,j)+ky3(k)*dt;  
            end
            for k=1:Nosc
                mod = Mod2(ax3(k),ay3(k));
                kx4(k) = FLdaOmegaCub(ax3(k),ay3(k),mod,lda(k),b(k),omega(k),a(k));
                kx4(k) = kx4(k)+sum(alpha(k,:).*ax3(:)');
                kx4(k) = kx4(k)-sum(alphaimag(k,:).*ay3(:)');
                ky4(k) = GLdaOmegaCub(ax3(k),ay3(k),mod,ldah(k),bh(k),omegah(k),ah(k));
                ky4(k) = ky4(k)+sum(beta(k,:).*ay3(:)');
                ky4(k) = ky4(k)+sum(betaimag(k,:).*ax3(:)');
            end
            for k=1:Nosc
                x(k,j+1) = x(k,j)+(kx1(k)+2*kx2(k)+2*kx3(k)+kx4(k))*dt/6;
                y(k,j+1) = y(k,j)+(ky1(k)+2*ky2(k)+2*ky3(k)+ky4(k))*dt/6;
            end
        end
    elseif NSETYPE == 1               
        for j=1:length(t)-1
            for k=1:Nosc
                if t(j) < ton(k)
                    x(k,j) = x(k,1);
                    y(k,j) = y(k,1);
                end
            end     
            if j == floor(tconn/dt)
                alpha = alphaconn;
                beta = betaconn;
                alphaimag = alphaimagconn;
                betaimag = betaimagconn;
            end
            for k=1:Nosc
                mod = Mod2(x(k,j),y(k,j));
                kx1(k) = FLdaOmegaCub(x(k,j),y(k,j),mod,lda(k),b(k),omega(k),a(k))+Fx(j);
                kx1(k) = kx1(k)+sum(alpha(k,:).*x(:,j)');
                kx1(k) = kx1(k)-sum(alphaimag(k,:).*y(:,j)');
                ky1(k) = GLdaOmegaCub(x(k,j),y(k,j),mod,ldah(k),bh(k),omegah(k),ah(k))+Fy(j); 
                ky1(k) = ky1(k)+sum(beta(k,:).*y(:,j)');
                ky1(k) = ky1(k)+sum(betaimag(k,:).*x(:,j)');
            end
            for k=1:Nosc
                ax1(k) = x(k,j)+kx1(k)*dt/2;    
                ay1(k) = y(k,j)+ky1(k)*dt/2;  
            end
            for k=1:Nosc
                mod = Mod2(ax1(k),ay1(k));
                kx2(k) = FLdaOmegaCub(ax1(k),ay1(k),mod,lda(k),b(k),omega(k),a(k))+(Fx(j)+Fx(j+1))/2;
                kx2(k) = kx2(k)+sum(alpha(k,:).*ax1(:)');
                kx2(k) = kx2(k)-sum(alphaimag(k,:).*ay1(:)');
                ky2(k) = GLdaOmegaCub(ax1(k),ay1(k),mod,ldah(k),bh(k),omegah(k),ah(k))+(Fy(j)+Fy(j+1))/2;
                ky2(k) = ky2(k)+sum(beta(k,:).*ay1(:)');
                ky2(k) = ky2(k)+sum(betaimag(k,:).*ax1(:)');
            end
            for k=1:Nosc
                ax2(k) = x(k,j)+kx2(k)*dt/2;    
                ay2(k) = y(k,j)+ky2(k)*dt/2;  
            end
            for k=1:Nosc
                mod = Mod2(ax2(k),ay2(k));
                kx3(k) = FLdaOmegaCub(ax2(k),ay2(k),mod,lda(k),b(k),omega(k),a(k))+(Fx(j)+Fx(j+1))/2;
                kx3(k) = kx3(k)+sum(alpha(k,:).*ax2(:)');
                kx3(k) = kx3(k)-sum(alphaimag(k,:).*ay2(:)');
                ky3(k) = GLdaOmegaCub(ax2(k),ay2(k),mod,ldah(k),bh(k),omegah(k),ah(k))+(Fy(j)+Fy(j+1))/2;
                ky3(k) = ky3(k)+sum(beta(k,:).*ay2(:)');
                ky3(k) = ky3(k)+sum(betaimag(k,:).*ax2(:)');
            end
            for k=1:Nosc
                ax3(k) = x(k,j)+kx3(k)*dt;    
                ay3(k) = y(k,j)+ky3(k)*dt;  
            end
            for k=1:Nosc
                mod = Mod2(ax3(k),ay3(k));
                kx4(k) = FLdaOmegaCub(ax3(k),ay3(k),mod,lda(k),b(k),omega(k),a(k))+Fx(j+1);
                kx4(k) = kx4(k)+sum(alpha(k,:).*ax3(:)');
                kx4(k) = kx3(k)-sum(alphaimag(k,:).*ay3(:)');
                ky4(k) = GLdaOmegaCub(ax3(k),ay3(k),mod,ldah(k),bh(k),omegah(k),ah(k))+Fy(j+1);
                ky4(k) = ky4(k)+sum(beta(k,:).*ay3(:)');
                ky4(k) = ky4(k)+sum(betaimag(k,:).*ax3(:)');
            end
            for k=1:Nosc
                x(k,j+1) = x(k,j)+(kx1(k)+2*kx2(k)+2*kx3(k)+kx4(k))*dt/6;
                y(k,j+1) = y(k,j)+(ky1(k)+2*ky2(k)+2*ky3(k)+ky4(k))*dt/6;
            end
        end
    elseif NSETYPE == 3        
        for j=1:length(t)-1
            for k=1:Nosc
                if t(j) < ton(k)
                    x(k,j) = x(k,1);
                    y(k,j) = y(k,1);
                end
            end     
            if j == floor(tconn/dt)
                alpha = alphaconn;
                beta = betaconn;
                alphaimag = alphaimagconn;
                betaimag = betaimagconn;
            end
            for k=1:Nosc
                mod = Mod2(x(k,j),y(k,j));
                kx1(k) = FLdaOmegaCub(x(k,j),y(k,j),mod,lda(k),b(k),omega(k),a(k))+Fx(j);
                kx1(k) = kx1(k)+sum(alpha(k,:).*x(:,j)');
                kx1(k) = kx1(k)-sum(alphaimag(k,:).*y(:,j)');
                ky1(k) = GLdaOmegaCub(x(k,j),y(k,j),mod,ldah(k),bh(k),omegah(k),ah(k))+Fy(j); 
                ky1(k) = ky1(k)+sum(beta(k,:).*y(:,j)');
                ky1(k) = ky1(k)+sum(betaimag(k,:).*x(:,j)');
            end
            for k=1:Nosc
                ax1(k) = x(k,j)+kx1(k)*dt/2;    
                ay1(k) = y(k,j)+ky1(k)*dt/2;  
            end
            for k=1:Nosc
                mod = Mod2(ax1(k),ay1(k));
                kx2(k) = FLdaOmegaCub(ax1(k),ay1(k),mod,lda(k),b(k),omega(k),a(k))+(Fx(j)+Fx(j+1))/2;
                kx2(k) = kx2(k)+sum(alpha(k,:).*ax1(:)');
                kx2(k) = kx2(k)-sum(alphaimag(k,:).*ay1(:)');
                ky2(k) = GLdaOmegaCub(ax1(k),ay1(k),mod,ldah(k),bh(k),omegah(k),ah(k))+(Fy(j)+Fy(j+1))/2;
                ky2(k) = ky2(k)+sum(beta(k,:).*ay1(:)');
                ky2(k) = ky2(k)+sum(betaimag(k,:).*ax1(:)');
            end
            for k=1:Nosc
                ax2(k) = x(k,j)+kx2(k)*dt/2;    
                ay2(k) = y(k,j)+ky2(k)*dt/2;  
            end
            for k=1:Nosc
                mod = Mod2(ax2(k),ay2(k));
                kx3(k) = FLdaOmegaCub(ax2(k),ay2(k),mod,lda(k),b(k),omega(k),a(k))+(Fx(j)+Fx(j+1))/2;
                kx3(k) = kx3(k)+sum(alpha(k,:).*ax2(:)');
                kx3(k) = kx3(k)-sum(alphaimag(k,:).*ay2(:)');
                ky3(k) = GLdaOmegaCub(ax2(k),ay2(k),mod,ldah(k),bh(k),omegah(k),ah(k))+(Fy(j)+Fy(j+1))/2;
                ky3(k) = ky3(k)+sum(beta(k,:).*ay2(:)');
                ky3(k) = ky3(k)+sum(betaimag(k,:).*ax2(:)');
            end
            for k=1:Nosc
                ax3(k) = x(k,j)+kx3(k)*dt;    
                ay3(k) = y(k,j)+ky3(k)*dt;  
            end
            for k=1:Nosc
                mod = Mod2(ax3(k),ay3(k));
                kx4(k) = FLdaOmegaCub(ax3(k),ay3(k),mod,lda(k),b(k),omega(k),a(k))+Fx(j+1);
                kx4(k) = kx4(k)+sum(alpha(k,:).*ax3(:)');
                kx4(k) = kx3(k)-sum(alphaimag(k,:).*ay3(:)');
                ky4(k) = GLdaOmegaCub(ax3(k),ay3(k),mod,ldah(k),bh(k),omegah(k),ah(k))+Fy(j+1);
                ky4(k) = ky4(k)+sum(beta(k,:).*ay3(:)');
                ky4(k) = ky3(k)+sum(betaimag(k,:).*ax3(:)');
            end
            for k=1:Nosc
                x(k,j+1) = x(k,j)+(kx1(k)+2*kx2(k)+2*kx3(k)+kx4(k))*dt/6;
                y(k,j+1) = y(k,j)+(ky1(k)+2*ky2(k)+2*ky3(k)+ky4(k))*dt/6;
            end
        end
    end


end
