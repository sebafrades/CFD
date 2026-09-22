% Se definen valores 

h=1;
N=8;
deltax=h/N;

% Relación L1/h

R1=1.25;

% Relación H/h

R2=1.75;

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

for i=1:(R1*N+1)
    A((R2-1)*(N/2)+1,i)=1;
end

for i=(R1*N+1):(R1*N+1)+(R3*N)
    A(1,i)=1;
end

for i=1:((R2-1)*(N/2)+1)
    A(i,(R1*N+1))=1;
end

for j=((R2-1)*(N/2)+1):(R2*N+1)-((R2-1)*(N/2)+1)
    A(j,1)=1-(j-((R2-1)*(N/2)+1))/N;
end

for j=1:R2*N+1
    A(j,(R1+R3)*N+1)=1-(j-1)/(R2*N);
end

% Se recorre la matriz 

for k=1:5000

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

% Se crea la matriz con los valores de velocidad V

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


quiver(X,Y,U,V)


streamline(X,Y,U,V,startx,starty)



% Se crea la matriz de eta (vorticidad)

B=zeros(round(R2*N+1),round((R1+R3)*N+1));

% Columna izquierda

for i=(R2-1)*N/2+2:(R2*N+1)-((R2-1)*(N/2)+1)
    B(i,1)=(V(i,2)-V(i,1))/deltax-(U(i-1,1)-U(i+1,1))/(2*deltax);
end

% Parte central

for i=(R2-1)*N/2+2:(R2*N+1)-((R2-1)*(N/2)+1) % fila 5 hasta 11
    for j=2:(R1+R3)*N % columna 2 hasta 18
        B(i,j)=(V(i,j+1)-V(i,j-1))/(2*deltax)-(U(i-1,j)-U(i+1,j))/(2*deltax);
    end
end

% Contorno izquierda, superior e inferior

for j=2:R1*N+1 % para columnas 2 a 11
    i=(R2-1)*N/2+1; % fila numero 4
    B(i,j)=(V(i,j+1)-V(i,j-1))/(2*deltax)-(U(i,j)-U(i+1,j))/deltax;
    i=(R2-1)*N/2+1+N; % fila numero 12
    B(i,j)=(V(i,j+1)-V(i,j-1))/(2*deltax)-(U(i-1,j)-U(i,j))/deltax;
end

for i=2:(R2*N) % filas 2 a 14
    for j=(R1*N+2):(R1+R3)*N % columnas 12 a 18
        B(1,j)=(V(1,j+1)-V(1,j-1))/(2*deltax)-(U(1,j)-U(2,j))/deltax;
        B((R2*N)+1,j)=(V((R2*N)+1,j+1)-V((R2*N)+1,j-1))/(2*deltax)-(U((R2*N),j)-U((R2*N)+1,j))/deltax;
        B(i,j)=(V(i,j+1)-V(i,j-1))/(2*deltax)-(U(i-1,j)-U(i+1,j))/(2*deltax);
    end
end

% para la columna derecha

for j=2:(R2*N) % desde 2 hasta 14
    B(j,(R1+R3)*N+1)=(V(j,(R1+R3)*N+1)-V(j,(R1+R3)*N))/(2*deltax)-(U(j-1,(R1+R3)*N+1)-U(j+1,(R1+R3)*N+1))/(2*deltax);
end

% Para los extremos
i=(R2-1)*N/2+1; % fila numero 4
B(i,1)=(V(i,2)-V(i,1))/deltax-(U(i,1)-U(i+1,1))/deltax;
i=(R2*N+1)-((R2-1)*(N/2)); %fila numero 12
B(i,1)=(V(i,2)-V(i,1))/deltax-(U(i-1,1)-U(i,1))/deltax;
j=(R1+R3)*N+1; % columna 19
B(1,j)=(V(1,j)-V(1,j-1))/deltax-(U(1,j)-U(2,j))/deltax;
i=(R2*N+1); %fila 15
B(i,j)=(V(i,j)-V(i,j-1))/deltax-(U(i-1,j)-U(i,j))/deltax;


