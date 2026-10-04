function [flag,Ps,Qs]=SepPlane(Ls,Bs,xb,yb)
nx=length(xb);ny=length(yb);
nbs=zeros(1,length(Bs));
for i=1:length(Bs)
    B=Bs{i}(xb,yb);nbs(i)=size(B,1);
end
nls=zeros(1,length(Ls));
for i=1:length(Ls)
    L=Ls{i}(xb,yb);nls(i)=size(L,1);
end
yalmip('clear')
Ps=cell(1,length(Ls));Qs=cell(1,length(Bs));
for i=1:length(Ps)
    Ps{i}=sdpvar(nls(i));
end
for i=1:length(Qs)
    Qs{i}=sdpvar(nbs(i));
end
zx=zeros(nx,1);zy=zeros(ny,1);
F=[];
for i=1:nx
    ex=zeros(nx,1);ex(i)=1;
    for j=1:ny        
        ey=zeros(ny,1);ey(j)=1;
        q=0;
        for k=1:length(Bs)            
            B=Bs{k}(ex,ey)-Bs{k}(ex,zy)-Bs{k}(zx,ey)+Bs{k}(zx,zy);
            q=q+trace(Qs{k}*B);
        end        
        F=[F,q==0];
    end    
end
w=0;
for i=1:length(Bs)
    F=[F,Qs{i}>=0];
    w=w+trace(Qs{i}*Bs{i}(xb,yb));
end
for i=1:length(Ls)
    F=[F,Ps{i}>=0];
    w=w+trace(Ps{i}*Ls{i}(xb,yb));
end
F=[F,w<=1e2];
op=sdpsettings;op.verbose=0;op.solver='SDPT3';
res=optimize(F,-w,op);          
for i=1:length(Ps)
    Ps{i}=value(Ps{i});
end
for i=1:length(Qs)
    Qs{i}=value(Qs{i});
end
flag=1;
zx=zeros(nx,1);zy=zeros(ny,1);
for i=1:nx
    ex=zeros(nx,1);ex(i)=1;
    for j=1:ny        
        ey=zeros(ny,1);ey(j)=1;
        q=0;
        for k=1:length(Bs)            
            B=Bs{k}(ex,ey)-Bs{k}(ex,zy)-Bs{k}(zx,ey)+Bs{k}(zx,zy);
            q=q+trace(Qs{k}*B);
        end
        if abs(q)>=1e-4
            flag=0;
        end
    end    
end
w=0;
for i=1:length(Bs)   
    if min(real(eig(Qs{i})))<=-1e-4
        flag=0;
    end
    w=w+trace(Qs{i}*Bs{i}(xb,yb));
end
for i=1:length(Ls) 
    if min(real(eig(Ps{i})))<=-1e-4
        flag=0;
    end
    w=w+trace(Ps{i}*Ls{i}(xb,yb));
end
if w<=1e-5
    flag=0;
end
end
