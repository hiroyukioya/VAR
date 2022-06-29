function  [INIT1,INIT2,INIT3,LIK,A,nR,NIS,Pshat]=TVVAR_KALMAN_FIS_R(rdata,order,T,INIT1,INIT2,INIT3)


data=detrend(rdata,'constant');

[tleng,nchan]=size(data); 

alpha=0.05;
beta=T;

% Prepare the Matrices for the State-Space model----------------------

% State transition matrix
F=eye(order*nchan^2);
% Driving matrix
G=eye(order*nchan^2); 
% Measurement noise covariance
pR=INIT3;
R=pR(:,:,ones(1,tleng));
% Process noise covariance
Q=eye(order*nchan^2)*beta; 
Q=Q(:,:,ones(1,tleng));
% Covariance matrix of the initial state vector
Pinit=INIT1;
% Covariance matrix of state vector estimate
Phat=Pinit(:,:,ones(1,tleng));
Pshat=Pinit(:,:,ones(1,tleng));
Pplus=Pinit(:,:,ones(1,tleng));
% Initial state vector
xinit=INIT2;
% State vector
xhat=zeros(nchan^2*order,tleng);
xshat=zeros(nchan^2*order,tleng);
xplus=zeros(nchan^2*order,tleng);
% normalized innovation squared
NIS=zeros(1,tleng);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                       Kalman Filter                                  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LIK=0;
for t=order+1:tleng
    if t==order+1
        Pminus=Pinit;
        xminus=xinit;
    end
    
  %  Measurement update.......................................
        % Observation matrix
        [H]=MakeObservationMatrix(data,t,nchan,order);
        % Innovation
        Inov=data(t,:)'-H*xminus;  
        % Covariance of Innovation
        V=H*Pminus*H'+R(:,:,t-1);
        % Kalman gain K
        invV=inv(V);
        Kgain=Pminus*H'*invV;    
        % State estimate
        xhat(:,t)=xminus+Kgain*Inov;
        % Error in the aposteriori estimate
        E=data(t,:)'-H*xhat(:,t);
        %  Covariance of error in the filtered estimate xhat
        Phat(:,:,t)=Pminus-Kgain*H*Pminus;
        %  Momentary covariance matrix of measurement noise (R2)
        R(:,:,t)=R(:,:,t-1)-alpha*(R(:,:,t-1)-E*E');
        %  Normalized Innovation Squared
        NIS(t)=Inov'*invV*Inov;
        %  Likelihood
        LIK=NIS(t)+log(abs(det(V)))+LIK;

   % Time update ......................................................  
        xplus(:,t)=F*xhat(:,t);
        Pplus(:,:,t)=F*Phat(:,:,t)*F'+G*Q(:,:,t)*G';
        xminus=xplus(:,t);
        Pminus=Pplus(:,:,t);
 end
 
 %********************************************%
 %%%%%        Kalman Fixed Interval Smoother           %%%%        
 %********************************************%
 
 for t=tleng:-1:order+1
     if t==tleng
           xsprior=xhat(:,t);
           Psprior=Phat(:,:,t);
     end
     %  Smoothing gain B
       B=Phat(:,:,t)*F'*inv(Pplus(:,:,t));
     %  Smoothed state vector
       xshat(:,t)=xhat(:,t)+B*(xsprior-xplus(:,t));
     %  Error covariance in state smoothing
       Pshat(:,:,t)=Phat(:,:,t)+B*(Psprior-Pplus(:,:,t))*B';
     %  Update xsprior and Psprior
       xsprior=xshat(:,t);
       Psprior=Pshat(:,:,t);
 end

B=xshat;
 
% Get initial values for next iteration
  INIT1=Pshat(:,:,100);
  INIT2=xshat(:,order+1);
  INIT3=mean(R(:,:,order+1:order+100),3);
  
% Get error covariance matrices  
  for t=order+1:tleng
        [H]=MakeObservationMatrix(data,t,nchan,order);
        Yhat(:,t)=H*B(:,t);
        I(:,t)=data(t,:)'-Yhat(:,t);  
        nR(:,:,t)=nR(:,:,t-1)-gamma*(nR(:,:,t-1)-I(:,t)*I(:,t)');
  end

     
LIK=-LIK;

 %  Rearrange the AR coefficients
       A=zeros(nchan,nchan,order,tleng);
       for t=1:tleng
             preA=reshape(B(:,t),order,nchan^2)';
             COEF=reshape(preA,nchan,nchan*order)'; 
             for ch=1:order
                   lm=1+nchan*(ch-1):nchan*ch;
                   A(:,:,ch,t)=COEF(lm,:);
             end
       end
  % Coefficients Matrices     
  A=permute(A,[1 2 4 3]);   
  
%------------------------------------------------------------------
function [H]=MakeObservationMatrix(data,m,nchan,order);
        for n=1:order
              ad(n)=m-n;
        end
        Temp=data(ad,:)';
        clear n
        for n=1:nchan
              preH(1+order*(n-1):order*n)=Temp(n,:);
        end
        I=eye(nchan);
        H=kron(I,preH);  % Kronecker product
