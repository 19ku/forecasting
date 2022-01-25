function x = adjout(y,thr,tflag)

%/* -- Adjust for outliers using fraction of IQR
%     
%      y = Data series
%      thr = threshold in multiples of IQR
%      tflag = 0  == replace with missing value 
%              1  == replace with maximum value
%              2  == replace with median value
%              3  == replace with local median (obs + or - 3 on each side)
%*

small = 1.0e-06;

missc=1e+32;
missv=missc.*ones(height(y),1);
%@ -- Compute IQR -- @
z=rmmissing(y,1);
z=sort(z,1);
zm=z(int16(height(z)/2)); 
iqr=z(int16(height(z)*3/4))-z(int16(height(z)/4)); %iqr=z(.75*height(z))-z(.25*height(z));

if iqr < small
    ya=standardizeMissing(0,0);
    x=ya;
    return
end

ya=abs(y-zm);

iya = gt(ya, (thr*iqr));
iyb = le(ya, (thr*iqr));
if tflag == 0
    x=(iyb .* y) + (iya .* missv);
    x=standardizeMissing(x,missc);
elseif tflag == 1
    isign = y > 0;
    jsign = -(y < 0);
    isign=isign+jsign;
    yt=(zm.*ones(rows(y),1)) + isign .* (thr .* ones(rows(y),1)); 
    x=(iyb .* y) + (iya .* yt);
elseif tflag == 2
    x=(iyb .* y) + (iya .* zm);
elseif tflag == 3
        %@ Compute rolling median @
    iwin=3;  %@ Window on either side @
    ymvec=standardizeMissing(zeros(rows(y),1),0);
    for i=1:height(y)
        j1=max([1,(i-iwin)]);
        j2=min([height(y),(i+iwin)]);
        x=rmmissing(y(j1:j2),1);
        x=sort(x,1);
        ymvec(i)=x(0.5*height(x));
    end
    x=(iyb .* y) + (iya .* ymvec);
    
end
end







