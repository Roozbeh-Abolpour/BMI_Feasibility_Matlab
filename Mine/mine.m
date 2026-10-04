function [X,Y]=mine(Ps)
n=length(Ps{1});r=length(Ps);
He=@(X) (X+X');
X=@(x) symreshape(x,n);
Y=@(y) reshape(y,n,n);
Bs=cell(1,r);
for i=1:r
    Bs{i}=@(x,y) Ps{i}+He(X(x)*Y(y));
end
Ls{1}=@(x,y) [-X(x)];
Ls{2}=@(x,y) [X(x)-eye(n)];
Ls{3}=@(x,y) [-eye(n) Y(y);Y(y)' -X(x)];
IP=InitPoint(n);
[xb,yb]=OuterLoop(Ls,Bs,IP);
X=X(xb);Y=Y(yb);
end
