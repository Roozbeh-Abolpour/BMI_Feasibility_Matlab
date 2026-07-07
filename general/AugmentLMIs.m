function [Ls,Bs]=AugmentLMIs(Ls,Bs,xb,yb,sign)
nx=length(xb);ny=length(yb);
for i=1:length(Ls)
    Ls{i}=@(x,y) Ls{i}(x,y(1:end-1));
end
for i=1:length(Bs)
    Bs{i}=@(x,y) Bs{i}(x,y(1:end-1));
end
Ls{end+1}=@(x,y) [-y(end)*eye(nx+ny) [x-xb;y(1:end-1)-yb];[x-xb;y(1:end-1)-yb]' -y(end)];
if sign==1
    Ls{end+1}=@(x,y) -1e4*(sum(x-xb)+sum(y(1:end-1)-yb))+y(end);
else
    Ls{end+1}=@(x,y) 1e4*(sum(x-xb)+sum(y(1:end-1)-yb))+y(end);
end
Ls{end+1}=@(x,y) -y(end)+1e-2;
end