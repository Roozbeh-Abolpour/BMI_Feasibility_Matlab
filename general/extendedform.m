function [Lss,Bss]=extendedform(Ls,Bs,nx,ny)
nl=length(Ls);nb=length(Bs);
Lss=cell(1,nl);Bss=cell(1,nb);
for k=1:nl
    Lss{k}=cell(1,1+nx+ny);
end
for k=1:nb
    Bss{k}=cell(1,1+nx+ny+nx*ny);
end
Ix=eye(nx);Iy=eye(ny);
for k=1:nl
    x=0*Ix(:,1);y=0*Iy(:,1);
    Lss{k}{1}=Ls{k}(x,y);
    for i=1:nx
        x=Ix(:,i);
        Lss{k}{1+i}=Ls{k}(x,y)-Lss{k}{1};
    end
    x=0*Ix(:,1);
    for i=1:ny
        y=Iy(:,i);
        Lss{k}{1+nx+i}=Ls{k}(x,y)-Lss{k}{1};
    end    
end
for k=1:nb
    x=0*Ix(:,1);y=0*Iy(:,1);
    Bss{k}{1}=Bs{k}(x,y);
    for i=1:nx
        x=Ix(:,i);
        Bss{k}{1+i}=Bs{k}(x,y)-Bss{k}{1};
    end
    x=0*Ix(:,1);
    for i=1:ny
        y=Iy(:,i);
        Bss{k}{1+nx+i}=Bs{k}(x,y)-Bss{k}{1};
    end    
    for i=1:nx
        x=Ix(:,i);
        for j=1:ny
            y=Iy(:,j);
            Bss{k}{1+nx+ny+(i-1)*ny+j}=Bs{k}(x,y)-...
                Bss{k}{1+i}-Bss{k}{1+nx+j}-Bss{k}{1};
        end
    end
end
end