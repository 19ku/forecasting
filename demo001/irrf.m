function [Impulse2] =irrf(ttt,k,reb,b,omega)

%irf    
C=zeros(k, ttt*k);
A=zeros(k, ttt*k);
C(1:k,1:k)=eye(3);
tb=reb';
tb=tb(:,2:end);
A(:,1:length(b)-1)=tb;

for i=1:ttt-1
    for j=1:i
        C(:,k*i+1:k*(i+1))=C(:,k*i+1:k*(i+1))+C(:,k*(i-j)+1:k*(i-j+1))*A(:,(j-1)*k+1:j*k);
    end
end

P=chol(omega)';   
Omega2=diag(diag(P));
B0=Omega2*inv(P);


irfx2=zeros(k,k*ttt);
for i=1:ttt
    irfx2(:,k*(i-1)+1:k*i) = C(:,k*(i-1)+1:k*i)*P;
end
    
Impulse2=zeros(ttt,k*k);
for j=1:ttt
    for i=1:k
        Impulse2(j,(i-1)*k+1:i*k)=irfx2(i,(j-1)*k+1:j*k);
    end
end
