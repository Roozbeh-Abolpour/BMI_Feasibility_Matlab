function x=syminv(X)
n=size(X,1);
x=zeros(n*(n+1)/2,1);it=1;
for i=1:n
    for j=i:n
        x(it)=X(i,j);it=it+1;
    end
end