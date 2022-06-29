function  [INIT2,INIT3,LIK,A,R,Inov,V,Q,Phat,xhat,NIS,AIC]=TVVAR_KALMAN_KR(rdata,order,T,INIT2,INIT3)


data=detrend(rdata,'constant');

[tleng,nchan]=size(data); 

alpha=0.03;
beta=T(1);

% Prepare the Matrices for the State-Space model----------------------

% State transition matrix
F=eye(order*nchan^2)*1;
% Driving matrix
G=eye(order*nchan^2);   
% Measurement noise covariance
pR=INIT3;
R=pR(:,:,ones(1,tleng));
% Process noise covariance
Q=eye(order*nchan^2)*beta; 
Q=Q(:,:,ones(1,tleng));
% Covariance matrix of the initial state vector
Pinit=eye(nchan^2*order)*10^-6;
% Initial state vector
xinit=INIT2;

H=zeros(nchan,nchan^2*order,tleng);
xhat=zeros(nchan^2*order,tleng);
Inov=zeros(nchan,tleng);
V=zeros(nchan,nchan,tleng);

% Make Observation Matrix H-------------------------------------
for m=order+1:tleng
        for n=1:order
              ad(n)=m-n;
        end
        Temp=data(ad,:)';
        clear n
        for n=1:nchan
              preH(1+order*(n-1):order*n)=Temp(n,:);
        end
        I=eye(nchan);
        H(:,:,m)=kron(I,preH);  % Kronecker product
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                  Run Kalman Filter                                %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LIK=0;
for t=order+1:tleng
    if t==order+1
        Pminus=Pinit;
        xminus=xinit;
    end
    
  %  Measurement update.......................................
        %  Innovation
        Inov(:,t)=data(t,:)'-H(:,:,t)*xminus;  
        %  Covariance matrix of the Innovation 
        V(:,:,t)=H(:,:,t)*Pminus*H(:,:,t)'+R(:,:,t-1);
        %  Kalman gain K
        Kgain=Pminus*H(:,:,t)'*inv(V(:,:,t));    
        %  State estimate
        xhat(:,t)=xminus+Kgain*Inov(:,t);
        %  Covariance of error in the filtered estimate xhat
        Phat=Pminus-Kgain*H(:,:,t)*Pminus;
        %  Momentary covariance matrix of measurement noise (R2)
        R(:,:,t)=R(:,:,t-1)-alpha*(R(:,:,t-1)-Inov(:,t)*Inov(:,t)');
        %  Normalized Innovation Squared
        NIS(t)=Inov(:,t)'*inv(V(:,:,t))*Inov(:,t);
        %  Likelihood
        LIK=NIS(t)+log(abs(det(V(:,:,t))))+LIK;

   % Time update ......................................................  
        xplus=F*xhat(:,t);
        Pplus=F*Phat*F'+G*Q(:,:,t)*G';
        xminus=xplus;
        Pminus=Pplus;
 end
 
 %  Rearrange the AR coefficients (= xshat)
       A=zeros(nchan,nchan,order,tleng);
       for t=1:tleng
             preA=reshape(xhat(:,t),order,nchan^2)';
             B=reshape(preA,nchan,nchan*order)'; 
             for ch=1:order
                   lm=1+nchan*(ch-1):nchan*ch;
                   A(:,:,ch,t)=B(lm,:);
             end
       end
  A=permute(A,[1 2 4 3]);           
  N=tleng-order;
  K=nchan^2*order*N;
  
  INIT2=xhat(:,end);
  INIT3=R(:,:,end);

%   AIC=2*LIK+log(size(rdata,1)-order)*order*nchan^2/(size(rdata,1)-order); 
%    AIC=2*LIK+2*log(tleng-order)*order*nchan^2*tleng/(tleng-order*nchan^2*tleng-1);
  %AIC=2*LIK+order*nchan^2*log(tleng-order);
%   AIC=2*LIK+2*K*N/(N-K-1);
  AIC=2*LIK+K*log(N);