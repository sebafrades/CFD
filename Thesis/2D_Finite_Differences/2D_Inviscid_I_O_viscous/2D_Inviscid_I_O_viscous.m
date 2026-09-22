% Se definen valores 

Q=1;
h=1;
N=8;
deltax=h/N;

% Relación L1/h

R1=1.25;

% Relación H/h

R2=1.75;
H=R2*h;

% Relación L2/h

R3=1;

% Deltax

deltax=h/N;

% Se crea la matriz

A=zeros(round(R2*N+1),round((R1+R3)*N+1));

% Se agregan las condiciones de contorno 

%(R2-1)*N/2+1
% R1*N+1
% R3*N

for j=((R2-1)*(N/2)+1):(R2*N+1)-((R2-1)*(N/2))
    x=((R2-1)*(N/2)+1);
    y=(R2*N+1)-((R2-1)*(N/2));
    a=-(4*(0.75*Q/h)*x*y)/(x^2+y^2-2*x*y);
    b=(4*(0.75*Q/h)*(x+y))/(x^2+y^2-2*x*y);
    c=-(4*(0.75*Q/h))/(x^2+y^2-2*x*y);
    d=-(a*y+(b*y^2)/2+(c*y^3)/3);
    A(j,1)=-(a*j+(b*j^2)/2+(c*j^3)/3+d);
end

for j=1:R2*N+1
    x=1;
    y=R2*N+1;
    a=-(4*(0.75*Q/H)*x*y)/(x^2+y^2-2*x*y);
    b=(4*(0.75*Q/H)*(x+y))/(x^2+y^2-2*x*y);
    c=-(4*(0.75*Q/H))/(x^2+y^2-2*x*y);
    d=-(a*y+(b*y^2)/2+(c*y^3)/3);
    A(j,(R1+R3)*N+1)=-(a*j+(b*j^2)/2+(c*j^3)/3+d);
end

for j=2:R1*N+1
    A(((R2-1)*(N/2)+1),j)=A(((R2-1)*(N/2)+1),1);
end

for i=1:((R2-1)*(N/2)+1)
    A(i,R1*N+1)=A(((R2-1)*(N/2)+1),1);
end

for j=R1*N+1:(R1+R3)*N
    A(1,j)=A(((R2-1)*(N/2)+1),1);
end

% Se recorre la matriz 

for k=1:250
    
    for i=(R2-1)*N/2+2:(R2*N+1)-((R2-1)*(N/2)+1)
        for j=2:(R1+R3)*N
            A(i,j)=(A(i+1,j)+A(i-1,j)+A(i,j-1)+A(i,j+1))/4;
        end
    end

    for i=2:(R2*N)
        for j=(R1*N+2):(R1+R3)*N
            A(i,j)=(A(i+1,j)+A(i-1,j)+A(i,j-1)+A(i,j+1))/4;
        end
    end
end

% Se crea la matriz con los valores de velocidad f

V=zeros(round(R2*N+1),round((R1+R3)*N+1));

for i=(R2-1)*N/2+2:(R2*N+1)-((R2-1)*(N/2)+1)
    for j=2:(R1+R3)*N
        V(i,j)=(A(i,j+1)-A(i,j-1))/(2*deltax);
    end
end

for i=2:(R2*N)
    for j=(R1*N+2):(R1+R3)*N
        V(i,j)=(A(i,j+1)-A(i,j-1))/(2*deltax);
    end
end

U=zeros(round(R2*N+1),round((R1+R3)*N+1));

for i=(R2-1)*N/2+2:(R2*N+1)-((R2-1)*(N/2)+1)
    for j=1:(R1+R3)*N
        U(i,j)=(A(i+1,j)-A(i-1,j))/(2*deltax);
    end
end

for j=1:R1*N+1
    U((R2-1)*N/2+1,j)=-(A((R2-1)*N/2+1,j)-A((R2-1)*N/2+2,j))/deltax;
    U((R2-1)*N/2+1+N,j)=(A((R2-1)*N/2+1+N,j)-A((R2-1)*N/2+N,j))/deltax;
end

for i=2:(R2*N)
    for j=(R1*N+2):(R1+R3)*N+1
        U(1,j)=(A(2,j)-A(1,j))/deltax;
        U((R2*N)+1,j)=(A((R2*N)+1,j)-A((R2*N),j))/deltax;
        U(i,j)=(A(i+1,j)-A(i-1,j))/(2*deltax);
    end
end


U=-U;

[X,Y]=meshgrid(0:deltax:(R1+R3)*h,0:deltax:R2*h);
starty=((R2-1)*N/2)*deltax:deltax:(((R2-1)*N/2+1)+N-2)*deltax;
%starty=transpose(starty);
%starty=flip(starty,1);
startx=zeros(size(starty));
streamline(X,Y,U,V,startx,starty)

quiver(X,Y,U,V)


streamline(X,Y,U,V,startx,starty)

