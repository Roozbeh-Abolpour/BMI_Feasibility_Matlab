function x=master(Lss,Bss,bounds,Uss,Vss,us,vs)
lbx=bounds.lbx;ubx=bounds.ubx;
lby=bounds.lby;uby=bounds.uby;
nx=length(lbx);ny=length(lby);
nl=length(Lss);nb=length(Bss);
yalmip('clear')
x=sdpvar(nx,1);
g=sdpvar(1,1);
F=[x>=lbx,x<=ubx,g>=-1e4];
for k=1:length(Uss)    
    Us=Uss{k};Vs=Vss{k};u=us{k};v=vs{k};
    h=u'*lby-v'*uby;
    for i=1:nl
        h=h+trace(Us{i}*Lss{i}{1});
        for j=1:nx
            h=h+trace(Us{i}*Lss{i}{1+j})*x(j);
        end
    end
    for i=1:nb
        h=h+trace(Vs{i}*Bss{i}{1});
        for j=1:nx
            h=h+trace(Vs{i}*Bss{i}{1+j})*x(j);
        end        
    end
    F=[F,g>=h];
end
op=sdpsettings;op.verbose=0;op.solver='mosek';
res=optimize(F,g,op);
if res.problem==0
    x=value(x);
else
    x=[];
end
end