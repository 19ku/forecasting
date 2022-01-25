clear all
clc

rng(3623541)

sname="demo_mc_ar4_i1";

load b_aic_i1.dat;
bar= b_aic_i1;
load lam1_aic_i1.dat;
lam1= lam1_aic_i1;
load lam2_aic_i1.dat;
lam2= lam2_aic_i1;
load stdr1_aic_i1.dat;
stdr1= stdr1_aic_i1;
load stdr2_aic_i1.dat;
stdr2= stdr2_aic_i1;
load smplvec_i1.dat;
smplvec= smplvec_i1;
%load eez.dat;


%  Sample Period @
fyds=1959;   % First Year of Data Set @
fmds=1;      % First Month of Data Set @
lyds=2002;   % Last Year of Data Set @
lmds=12;     % Last Month of Data Set @
fyest=1960;  % First Year to For Estimation @
fmest=1;     % First Month For Estimation @
fyforc=1979; % First Year to For Forecasting @
fmforc=1;    % First Month For Forecasting @
lyforc=2002; % Last Year to For Forecasting @
lmforc=11;   % Last Month For Forecasting @
bdy=1983;    % Break Date Yearr @
bdm=1;       % Break Date Month @

% Other Parameters @
nrep=10;             % Actual Number of replications @
nit=100;            % Number of initial periods for startup @
nphvec=[3; 6; 12; 24];   % Number of periods ahead for forecast, >1 @
[rnphvec, cnphvec] = size(nphvec);
nestmin=120;        % Minimum number of Observations for estimation @
narest=4;          % Number of lags for estimation @

nphmax=max(nphvec);

thr=6;

ns=width(bar);
nar=height(bar)-1;

const=bar(1,:);
arcoef=bar(2:height(bar),:);
rev_arcoef=flip(arcoef);

%@ -- Check for Stationarity -- @
%i=1;
%while i <= cols(bar)
% b=arcoef(:,i);
% if maxc(abs(polyroot(1|-b))) > 1.0;
%  "Unstable Model for series";;i;
%  "Processing Stops";
%  break;
% i=i+1; 
%end
 %@ -- Sample Period Pointers -- @
nds =   12*(lyds-fyds)  + (lmds-fmds)+1;  %@ Data Set Sample Size @
ifest = 12*(fyest-fyds) + (fmest-fmds)+1; %@ Index of First Estimation Period @
ifforc = 12*(fyforc-fyds) + (fmforc-fmds)+1; %@ Index of First Forecast Period@
ilforc = 12*(lyforc-fyds) + (lmforc-fmds)+1; %@ Index of Last Forecast Period @
ifbd =  12*(bdy-fyds) + (bdm-fmds)+1;     %@ Index of First Forecast Period @
nds1 = ifbd-1;
nds2=nds-nds1;
nt1 = nds1 + nit;
nt2 = nds2;
nt=nt1+nt2;

 %@ -- Construct Calendar Sequences for Plotting, etc -- @ 
 calds=zeros(nds,2);
 calds(1,1)=fyds;
 calds(1,2)=fmds;
 yr=fyds; mt=fmds;
 i=2;
 while i <= nds
    mt=mt+1;
    if mt > 12
        mt=1; yr=yr+1;
    end
    calds(i,1)=yr; calds(i,2)=mt;
    i=i+1; 
 end
 calest=calds(ifest:height(calds),:);
 tds=calds(:,1)+(calds(:,2)-ones(height(calds),1))/12;
 test=calest(:,1)+(calest(:,2)-ones(height(calest),1))/12;
 
% @ -- Update Sample Pointers for Forecast Period -- @
ifforc=ifforc-ifest+1;
ilforc=ilforc-ifest+1;

%@ Construct Matrix with 0's and missing values consistent with smplvec @
trnd=linspace(1,1*nds,nds)';
tmat=trnd.*ones(1,ns);
tmat=(tmat >= smplvec(1,:)).*(tmat <= smplvec(2,:));
tmat=standardizeMissing(tmat,0);
tmat=ones(height(tmat),width(tmat))-tmat;

irep=1;
while irep<=nrep
    rsltvec=standardizeMissing(zeros(ns, height(nphvec)),0);
    
% @ -- Generate errors -- @
    fac1=randn(nt1, width(lam1));
    uniq1=randn(nt1,ns).*stdr1';
    e1=(fac1*lam1')+uniq1;

    fac2=randn(nt2, width(lam2));
    uniq2=randn(nt2,ns).*stdr2';
    e2=(fac2*lam2')+uniq2;
    
    e=[e1;e2];
    
% @ -- Generate Data -- @
    ydata=zeros(nt,ns);
    for t=nar+1:nt
        ydata(t,:)=const;
        ydata(t,:)=const+(sum(ydata(t-nar:t-1).*rev_arcoef))'+e(t,:);
    end
    ydata=ydata(nit+1:nt,1);
%@ Convert to Missing Values as determined by smplvec @
    ydata=ydata+tmat;

    is=1;
    while is <=width(ydata)
        yforrslt=standardizeMissing(zeros(height(test), 1),0);
        yforrsp=standardizeMissing(zeros(nphmax, 1),0);
        
        inph=1;
        while inph<=height(nphvec)
            syfor=strcat('pow', num2str(inph,"%02d"),'.txt');
            %fileID = fopen(syfor,'w');
            writematrix(yforrslt,syfor); 
            syfor=strcat('prj', num2str(inph,"%02d"),'.txt');
            %fileID = fopen(syfor,'w');
            writematrix(yforrslt,syfor);             
            inph=inph+1;
        end
        tcode=smplvec(3,is);
        y=ydata(:,is);
        %@ -- Adjust for outliers -- to mimic steps in empirical programs @
        ya=adjout(y,thr,0);
        y=ya; 
        
  %@ ----- Construct power up forecast ----- @  
   %@ -- Construct Matrix of Own lags to be used in forecasting -- @
        yreg=zeros(height(y),nar);
        yreg=standardizeMissing(yreg,0);
        i=1; 
        while i <= nar %@ yreg contains lags @
            yreg(i+1:height(yreg),i)=y(1:height(y)-i);
            i=i+1; 
        end
        %@ Set up Matrix of Control Variables @
        cvar=ones(height(yreg),1); %@ Constant Term @


        %@ -- Save only Period relevant for estimation/forecasting -- @
        yfor=y(ifest:height(y));
        yreg=yreg(ifest:height(yreg),:);
        cvar=cvar(ifest:height(cvar),:);
        
        t=ifforc;
        while t <= ilforc
            yfort=yfor(1:t);   %@ Dependent Variable in Regressions @ 
            yregt=yreg(1:t,:);
            cvart=cvar(1:t,:);
            
        %@ -- Eliminate any missing values -- @
            tmp=rmmissing([yfort, cvart, yregt],1);
            if height(tmp) < nestmin 
                %goto skipt1; 
                continue
            end
            yfort=tmp(:,1);
            cvart=tmp(:,2);
            yregt=tmp(:,3:width(tmp));
            
            %@ Run AutoRegression @
            xt=cvart;
            if narest > 0
                xt=[xt, yregt(:,1:narest)]; 
            end
            b=(inv(xt'*xt))*(xt'*yfort);
            dd=nar-narest;
            if dd > 0
            	%@ pad with zeros @
                b=[b;zeros(dd,1)];
            end
            
        %@ Compute Dynamic Forecasts of y @
            if nar > 0
            %@ Construct dynamic forecast @ 
                if nar > 1
                    id = eye(nar);
                    idr = id(1:nar-1,:);
                    compm = [b(2:height(b),:)' ; idr]; %@ companion form matrix @
                    cvarm = [b(1,:) ; idr(:,nar:nar)]; %@ companion form constant (1 0 0...0)'@
                else %@ nar=1 @
                    compm = b(2:rows(b),:); 
                    cvarm = b(1:1,:); 
                end
                %@ 1 to nphmax step ahead forecasts @
                hh=1;
                while hh<=nphmax
                    if hh==1 
                        compmp=compm;
                        cvarmp=cvarm;
                    else
                        cvarmp=cvarm+compm*cvarmp;
                        compmp=compm*compmp;
                    end
                    yforrsp(hh)= [cvar(t+1,:),yreg(t+1,1:nar)]*[cvarmp(1,1),compmp(1,:)]'; 
                    %@ hh step ahead forecast @
                hh=hh+1;
                end
            end
            if nar == 0
                yforrsp=b(1)*ones(nphmax,1);
            end
            
            %@ -- Construct Forecasts over Different Horizons -- @
            inph=1;
            while inph <= height(nphvec)
                nph=nphvec(inph);
                syfor=strcat('pow', num2str(inph,"%02d"),'.txt');
                yforrslt = readmatrix(syfor);
                
                %@ -- AR Forecast, power (aggregate forecasts to match proj definition) -- @
                if tcode==1 || tcode==4     %@ I(0) @
                    yforrslt(t,1)=yforrsp(nph);
                elseif tcode==2 || tcode==5 %@ I(1) @
                    yforrslt(t,1)=sum(yforrsp(1:nph));
                elseif tcode==3 || tcode==6 %@ I(2) @
                    yforrslt(t,1)=sum(cumsum(yforrsp(1:nph)));
                end 
                
                %@ -- temporarily save forecasts -- @
                syfor = strcat('pow', num2str(inph,"%02d"),'.txt');
                writematrix(yforrslt,syfor);        
                inph=inph+1; 
            end
      
            %skipt1;
    
            t=t+1; 
        end
        
        %@ ----- Construct direct forecasts ----- @  
  
        inph=1;
        while inph <= height(nphvec)
        nph=nphvec(inph);
        syfor = strcat('prj', num2str(inph,"%02d"),'.txt');
        yforrslt = readmatrix(syfor);
        yfor=yfcsta(y,tcode,nph); %@ -- Construct Forecast variable --
                                  %Note: yfor(t) is variable being forecast
                                  %at time t -- e.g., 
                                  %yfor(t)=x(t+nph)-x(t)
                                  %for tcode 2 (first differences) @
                                  
        %@ -- Construct Matrix of Own lags to be used in forecasting -- @
        yreg=zeros(height(y),nar);
        yreg=standardizeMissing(yreg,0);
        i=1;
        while i <= nar
            yreg(i:height(yreg),i)=y(1:height(y)+1-i);
        i=i+1;
        end
        %@ Set up Matrix of Control Variables @
        cvar=ones(height(yreg),1); %@ Constant Term @
        
        %@ -- Save only Period relevant for estimation/forecasting -- @
        yfor=yfor(ifest:height(yfor));
        yreg=yreg(ifest:height(yreg),:);
   
        t=ifforc;
        while t <= ilforc
            yfort=yfor(1:t-nph);   %@ Dependent Variable in Regressions @ 
            %fprintf('#%d\n',t);
            
            
            %@ AR Forecast @
            %@ Regressors @
            yregt=yreg(1:t-nph,:);
            cvart=cvar(1:t-nph,:);
      
            %@ -- Eliminate any missing values -- @
            tmp=rmmissing([yfort,cvart,yregt]);
            if height(tmp) < nestmin
              %goto skipt; 
	      t=t+1;
              continue
            end
            yfort=tmp(:,1);
            cvart=tmp(:,2);
            yregt=tmp(:,3:width(tmp));
            
           %@ Run AutoRegression @
            xt=cvart;
            if narest > 0
                xt=[xt, yregt(:,1:narest)]; 
            end
            b=(inv(xt'*xt))*(xt'*yfort);
            dd=nar-narest;
            if dd > 0
            %@ pad with zeros @
                b=[b;zeros(dd,1)];
            end
    
          %@ Construct AR Forecast @
            zt=cvar(t,:);
            if nar > 0
                zt=[zt, yreg(t,1:nar)]; 
            end
            yforrslt(t,1)=zt*b;  %@ -- AR Forecast -- @

        %skipt;
        t=t+1; 
        end
        
        %@ -- temporarily save forecasts -- @
        syfor = strcat('prj', num2str(inph,"%02d"),'.txt');
        writematrix(yforrslt,syfor);     
        inph=inph+1; 
        end
        
        
    %@ -- Compute RMSE for different Forecasting Horizons -- @
  
    inph=1; 
    while inph <= height(nphvec)
    nph=nphvec(inph);
    yfor=yfcsta(y,tcode,nph);  %@ Actual @
    yfor=yfor(ifest:height(yfor));
    syfor = strcat('prj', num2str(inph,"%02d"),'.txt');
    yforprj = readmatrix(syfor);
    syfor = strcat('pow', num2str(inph,"%02d"),'.txt');
    yforpow = readmatrix(syfor);  
    eprj=yfor-yforprj;
    epow=yfor-yforpow;
    ee=rmmissing([epow, eprj]);
    if height(ee) > 36
        e2=mean(ee.^2);
        relmse=e2(2)/e2(1);
        rsltvec(is,inph)=relmse; 
    end 
    inph=inph+1; 
    end
    is=is+1; 
    end  
    
    if irep == 1
        res_all=reshape(rsltvec,[],1);
    else
        res_all=[res_all, reshape(rsltvec,[],1)];
    end
    syfor = strcat(sname,'.txt');
    writematrix(res_all,syfor);     
    fprintf('Just finished iteration #%d\n', irep);
    irep=irep+1;
end