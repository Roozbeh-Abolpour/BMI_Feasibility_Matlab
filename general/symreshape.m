function X=symreshape(x,n)
it=1;
for i=1:n
    for j=i:n
        X(i,j)=x(it);it=it+1;
        X(j,i)=X(i,j);
    end
end
end