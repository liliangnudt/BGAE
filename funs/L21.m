function [E]=L21(Q,lambda)

E=zeros(size(Q));
for i=1:size(Q,2)
    temp=norm(Q(:,i),2);
    if temp<lambda
        E(:,i)=0;
    else
        E(:,i)=(temp-lambda)/temp*Q(:,i);
    end
end
