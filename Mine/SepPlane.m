function [Ps,Qs]=SepPlane(Ls,Bs,x,y,kappa)
nx=length(x);ny=length(y);
nl=length(Ls);nb=length(Bs);
yalmip('clear')
Ps=cell(1,nl);Qs=cell(1,nb);
F=[];
for i=1:nl
    m=length(Ls{i}(x,y));
    Ps{i}=sdpvar(m,m,'symmetric');
    F=[F,Ps{i}>=0];
end
for i=1:nb
    m=length(Bs{i}(x,y));
    Qs{i}=sdpvar(m,m,'symmetric');
    F=[F,Qs{i}>=0];
end
Ix=eye(nx);Iy=eye(ny);
for i=1:nx
    for j=1:ny
        R=0;
        for k=1:nb
            B=Bs{k}(Ix(:,i),Iy(:,j))-Bs{k}(Ix(:,i),0*Iy(:,j))-...
                Bs{k}(0*Ix(:,i),Iy(:,j))+Bs{k}(0*Ix(:,i),0*Iy(:,j));
            R=R+Qs{k}*B;
        end
        F=[F,trace(R)==0];
    end
end
delta=0;
for k=1:nl
    L=Ls{k}(x,y);
    delta=delta+trace(Ps{k}*L);
end
for k=1:nb
    B=Bs{k}(x,y);
    delta=delta+trace(Qs{k}*B);
end
F=[F,delta<=kappa];
op=sdpsettings;op.verbose=0;op.solver='mosek';
res=optimize(F,-delta,op);
if res.problem~=0
    Ps={};Qs={};
else
    for i=1:nl
        Ps{i}=value(Ps{i});
    end
    for i=1:nb
        Qs{i}=value(Qs{i});
    end
end
end