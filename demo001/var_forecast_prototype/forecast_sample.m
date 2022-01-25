clear all
clc

%load omochi.dat;

%omochi=xlsread('seminar_var_ex1_data.xlsx');  % 77x3 matrix
%omochi=xlsread('varseminar_data.xlsx');
load quart_last.dat;
omochi=quart_last;



lag=3;
[T, k]=size(omochi);



X=omochi(lag:T-1,:);
X0=omochi(1,:);
for i=1:lag-1
    X=[X omochi(lag-i:T-(i+1),:)];
    X0=[X0 omochi(i+1,:)];
end

x=omochi(lag+1:end,:);
X=[ones(T-lag,1), X];

b=inv(X'*X)*X'*x;

e=x-X*b;
omega=e'*e/(T-lag-lag*k-1);

forecast_horizon=12;
yhat_store=zeros(forecast_horizon,k);

capy=X(end,:)';
for i=1:forecast_horizon
    yhat=b'*capy;
    yhat_store(i,:)=yhat';
    capy=[1; yhat; capy(2:end-k,:)];
end    

nar= length(b)-1;
narr=(lag-1)*k;

id = kron(eye(lag),eye(k));
idr = id(1:narr,:);
compm = [b(2:height(b),:)' ;idr];
cvarm = [b(1,:)';idr(:,end:end)];

yforrsp=zeros(forecast_horizon,k);

for hh=1:forecast_horizon
    if hh==1
        compmp=compm;
        cvarmp=cvarm;
    else
        cvarmp=cvarm+compm*cvarmp;
        compmp=compm*compmp;
    end
    yforrsp(hh,:)=X(end,:)*[cvarmp(1:k,:), compmp(1:k,:)]';

end



limit_T=100;

for corp=T-limit_T:T-forecast_horizon


X=omochi(lag:corp-horizon_forecast,:);
X0=omochi(1,:);
for i=1:lag-1
    X=[X omochi(lag-i:T-(i+1),:)];
    X0=[X0 omochi(i+1,:)];
end

x=omochi(lag+1:corp+forecast_horizon,:);
X=[ones(T-lag,1), X];

Zt=X;


b=inv(X'*X)*X'*x;


yhat=b'*capy=X(end,:)';

capy=X(end,:)';
for i=1:forecast_horizon
    yhat=b'*capy;
    yhat_store(i,:)=yhat';
    capy=[1; yhat; capy(2:end-k,:)];
end    

end