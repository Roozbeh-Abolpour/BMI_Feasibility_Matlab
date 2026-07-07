function [x,y,W,g]=bablb(Lss,Bss,lbx,ubx,lby,uby)
nx=length(lbx);ny=length(lby);
nl=length(Lss);nb=length(Bss);
yalmip('clear')
x=sdpvar(nx,1);y=sdpvar(ny,1);
g=sdpvar(1,1);W=sdpvar(nx,ny,'full');
F=[x>=lbx,x<=ubx,y>=lby,y<=uby];
for i=1:nx
    for j=1:ny
        F=[F,lbx(i)*y(j)-W(i,j)-lbx(i)*lby(j)+lby(j)*x(i)<=0];
        F=[F,lbx(i)*uby(j)-x(i)*uby(j)-lbx(i)*y(j)+W(i,j)<=0];
        F=[F,x(i)*uby(j)-ubx(i)*uby(j)-W(i,j)+ubx(i)*y(j)<=0];
        F=[F,W(i,j)-ubx(i)*y(j)-x(i)*lby(j)+ubx(i)*lby(j)<=0];
    end
end
for k=1:nl
    L=Lss{k}{1};
    for i=1:nx
        L=L+Lss{k}{1+i}*x(i);
    end
    for j=1:ny
        L=L+Lss{k}{1+nx+j}*y(j);
    end    
    F=[F,L<=g*eye(length(L))];
end
for k=1:nb
    B=Bss{k}{1};
    for i=1:nx
        B=B+Bss{k}{1+i}*x(i);
    end
    for j=1:ny
        B=B+Bss{k}{1+nx+j}*y(j);
    end    
    for i=1:nx
        for j=1:ny
            B=B+Bss{k}{1+nx+ny+(i-1)*ny+j}*W(i,j);
        end
    end
    F=[F,B<=g*eye(length(B))];
end
op=sdpsettings;op.solver='mosek';op.verbose=0;
sol=optimize(F,g,op);
if sol.problem==0
    x=value(x);
    y=value(y);
    W=value(W);
    g=value(g);
else
    x=[];y=[];W=[];g=[];
end
end