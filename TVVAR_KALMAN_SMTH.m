function  [INIT1,INIT2,INIT3,LIK,A,R,V,VV,nR]=TVVAR_KALMAN_SMTH(rdata,order,T,INIT1,INIT2,INIT3,Fs)

% INIT2=" initial state vector "
% INIT3=" initial measurement covariance "


data=detrend(rdata,'constant');

[tleng,nchan]=size(data); 

alpha=0.03;
beta=T(1);

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
Pplus=Pinit(:,:,ones(1,tleng));
% Initial state vector
xinit=INIT2;
%state vector
xhat=zeros(nchan^2*order,tleng);
xshat=zeros(nchan^2*order,tleng);
xplus=zeros(nchan^2*order,tleng);
% Innovation cov.
Inov=zeros(nchan,tleng);
V=zeros(nchan,nchan,tleng);
VV=zeros(nchan,nchan,tleng);

% Make Observation Matrix H-------------------------------------
% H=zeros(nchan,nchan^2*order,tleng);
% for m=order+1:tleng
%         for n=1:order
%               ad(n)=m-n;
%         end
%         Temp=data(ad,:)';
%         clear n
%         for n=1:nchan
%               preH(1+order*(n-1):order*n)=Temp(n,:);
%         end
%         I=eye(nchan);
%         H(:,:,m)=kron(I,preH);  % Kronecker product
%  end
 
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
        %Observation matrix
        [H]=MakeObservationMatrix(data,t,nchan,order);
        %  Innovation
        Inov(:,t)=data(t,:)'-H*xminus;  
        %  Covariance matrix of the Innovation 
        V(:,:,t)=H*Pminus*H'+R(:,:,t-1);
        %  Kalman gain K
        Kgain=Pminus*H'*inv(V(:,:,t));    
        %  State estimate
        xhat(:,t)=xminus+Kgain*Inov(:,t);
        %  Covariance of error in the filtered estimate xhat
        Phat(:,:,t)=Pminus-Kgain*H*Pminus;
        % error
        ER=data(t,:)'-H*xhat(:,t);
        %  Momentary covariance matrix of measurement noise (R2)
        R(:,:,t)=R(:,:,t-1)-alpha*(R(:,:,t-1)-Inov(:,t)*Inov(:,t)');
       %  Normalized Innovation Squared
        NIS(t)=Inov(:,t)'*inv(V(:,:,t))*Inov(:,t);
        %  Likelihood
        LIK=NIS(t)+log(abs(det(V(:,:,t))))+LIK;

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
     %  Error covariance
       Pshat=Phat(:,:,t)+B*(Psprior-Pplus(:,:,t))*B';
     %  Update xsprior and Psprior
       xsprior=xshat(:,t);
       Psprior=Pshat;
 end

clear B;
B=xshat;

% Smoothing Coefficients
%    [B]=smothe(xshat,Fs);
%  for t=2:tleng
%        B(:,t)=B(:,t-1)-0.03*(B(:,t-1)-xshat(:,t));
%  end

 
  % Get initial values for next iteration
  INIT1=Pshat;
  INIT2=xshat(:,order+1);
  nR=R;
  for t=order+1:tleng
        [H]=MakeObservationMatrix(data,t,nchan,order);
        I(:,t)=data(t,:)'-H*B(:,t);  
        nR(:,:,t)=nR(:,:,t-1)-0.03*(nR(:,:,t-1)-I(:,t)*I(:,t)');
  end
%   INIT3=mean(nR(:,:,order+1:order+10),3);
    INIT3=mean(nR(:,:,end-10:end),3);
      
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
%------------------------------------------------------------------
function [G]=smothe(S,Fs);
            G=WindowedLowPassFilter(S',Fs,50,1);
            G=G';

function [new]=WindowedLowPassFilter(sig,Fs,N,Fc1);
flag = 'scale';  % Sampling Flag
win = hamming(N+1);
b  = fir1(N, [Fc1]/(Fs/2), 'low', win, flag);
Hd = dfilt.dffir(b);
win = hamming(2*N+1);
z  = fir1(2*N, [Fc1]/(Fs/2), 'low', win, flag);
[new]=filtfilt(b,1,sig);

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
