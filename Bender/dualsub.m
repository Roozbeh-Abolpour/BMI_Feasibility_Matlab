function [Us,Vs,u,v]=dualsub(Lss,Bss,bounds,x)
lbx=bounds.lbx;ubx=bounds.ubx;
lby=bounds.lby;uby=bounds.uby;
nx=length(lbx);ny=length(lby);
nl=length(Lss);nb=length(Bss);
yalmip('clear')
Us=cell(1,nl);Vs=cell(1,nb);
u=sdpvar(ny,1);v=sdpvar(ny,1);
F=[u>=0,v>=0];
h=0;
for i=1:nl
    m=size(Lss{i}{1},2);
    Us{i}=sdpvar(m,m,'symmetric');
    h=h+trace(Us{i});
    F=[F,Us{i}>=0];
end
for i=1:nb
    m=size(Bss{i}{1},2);
    Vs{i}=sdpvar(m,m,'symmetric');
    h=h+trace(Vs{i});
    F=[F,Vs{i}>=0];
end
F=[F,h==1];
for k=1:ny
    h=0;
    for i=1:nl
        h=h+trace(Us{i}*Lss{i}{1+nx+k});
    end
    for i=1:nb
        h=h+trace(Vs{i}*Bss{i}{1+nx+k});
        for j=1:nx
            h=h+trace(Vs{i}*Bss{i}{1+nx+ny+(j-1)*ny+k})*x(j);
        end
    end
    F=[F,h-u(k)+v(k)==0];
end
g=u'*lby-v'*uby;
for i=1:nl
    g=g+trace(Us{i}*Lss{i}{1});
    for j=1:nx
        g=g+trace(Us{i}*Lss{i}{1+j})*x(j);
    end
end
for i=1:nb
    g=g+trace(Vs{i}*Bss{i}{1});
    for j=1:nx
        g=g+trace(Vs{i}*Bss{i}{1+j})*x(j);
    end
end
ops=sdpsettings;ops.solver='mosek';ops.verbose=0;
sol=optimize(F,-g,ops);
if sol.problem==0
    u=value(u);v=value(v);
    for i=1:nl
        Us{i}=value(Us{i});        
    end
    for i=1:nb
        Vs{i}=value(Vs{i});        
    end    
else
    Us={};Vs={};u=[];v=[];
end
end