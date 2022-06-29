function  [LIK,R,Inov,V,Q,Phat,xhat]=TVVAR_KALMAN_LIK(T)
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

global ORD TNUM FNEW FIRST SECOND INIT1 INIT2 INIT3

rdata=FNEW(:,[FIRST SECOND],TNUM);
order=ORD;

rdata=detrend(rdata,'constant');
data=rdata;
%  [data]=varequalize_r(rdata);
[tleng,nchan]=size(data); 

alpha=0.05;
beta=T(1);

 [coef,INIT3]=VAR_EST_QR(data(1:120,:),order);
 [INIT2]=rearrange(coef,nchan,order);
%   INIT2=zeros(nchan^2*order,1);


 
% Prepare the Matrices for the State-Space model----------------------

% State transition matrix
F=eye(order*nchan^2)*1;
% Driving matrix
G=eye(order*nchan^2);   
% Measurement noise covariance
pR=INIT3;
R=pR(:,:,ones(1,tleng));
% Process (State) noise covariance
Q=eye(order*nchan^2)*beta; 
Q=Q(:,:,ones(1,tleng));
% Covariance matrix of the error in the initial state vector
Pinit=eye(nchan^2*order)*10^-6;
% Initial state vector
xinit=INIT2;

H=zeros(nchan,nchan^2*order,tleng);
xhat=zeros(nchan^2*order,tleng);
Inov=zeros(nchan,tleng);
V=zeros(nchan,nchan,tleng);
xplus=zeros(nchan^2*order,tleng);

% Prepare Observation Matrices H--------------------------------------
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
        xplus(:,t)=F*xhat(:,t);
        Pminus=F*Phat*F'+G*Q(:,:,t)*G';
        xminus=xplus(:,t);
 end
   
% %  Calculate log likelihood...........................................
%         LIK=0;
%         st=order;
%         en=tleng;
%   for n=st+1:en
%         preINOa=NIS(n);
%         preSUV=log((det(V(:,:,n))));
%         preLIK=preINOa+preSUV;
%         LIK=preLIK+LIK;
%   end
  
   