function [x,g]=babubx(y,Ls,Bs,lbx,ubx)
nx=length(lbx);
yalmip('clear')
x=sdpvar(nx,1);
g=sdpvar(1,1);
F=[x>=lbx,x<=ubx];
for i=1:length(Ls)
    L=Ls{i}(x,y);
    F=[F,L<=g*eye(length(L))];
end
for i=1:length(Bs)
    B=Bs{i}(x,y);
    F=[F,B<=g*eye(length(B))];
end
ops=sdpsettings;ops.solver='mosek';ops.verbose=0;
sol=optimize(F,g,ops);
if sol.problem==0
    g=value(g);x=value(x);
else
    x=[];g=inf;   
end
end