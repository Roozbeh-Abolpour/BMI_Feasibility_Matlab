function g=ErrBound(Ls,Bs,x,y)
g=-inf;
for i=1:length(Ls)
    L=Ls{i}(x,y);g=max(real(eig(L)));
end
for i=1:length(Bs)
    B=Bs{i}(x,y);g=max(real(eig(B)));
end
end