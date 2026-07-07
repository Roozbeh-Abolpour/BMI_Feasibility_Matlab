function ef=checksol(Ls,Bs,x,y)
if isempty(x)
    ef=-1;
    return
end
ef=1;
for i=1:length(Ls)
    L=Ls{i}(x,y);
    if max(real((eig(L))))>1e-4
        ef=0;
        return
    end
end
for i=1:length(Bs)
    B=Bs{i}(x,y);
    if max(real((eig(B))))>1e-4
        ef=0;
        return
    end
end
end