%% Función vorticidad - líneas de corriente, expansión abrupta de un flujo

%%

% Se definen valores 

H=2;
N=15; % número de nodos
iteraciones=2;
deltax=(H/2)/(N-1); 
Re=400;
MaxIt=300000;
maxe = 1e-7;

% Relación H/h

R2=2; 
h=H/R2;

% Relación L1/h

R1=1;

% Relación L2/h

R3=20;

% Se crea la matriz de líneas de corriente

cfinal=(R1/R2+R3/R2)*2*(N-1)+1; % número de columnas %(R1+R3)*(N-1)+1

%% Se crea la matriz A, líneas de corriente, con las condiciones de contorno y se resuelve el caso potencial

Aneu=zeros(N,cfinal);

% Se agregan las condiciones de contorno 

% Columna izquierda, régimen Poiseuille

y=linspace(H,H/2,N);

fh=ceil(2*(N-1)*(0.5-1/(2*R2))+0.5); % fila del valor de h %((1-1/R2)*(N-1)+1)

for j=fh:N % para la fila 4 hasta 8

    y4=y(fh);

    a=2/(H^3-8*y4^3+12*H*y4^2-6*H^2*y4);
    b=(-3*H)/(H^3-8*y4^3+12*H*y4^2-6*H^2*y4);
    c=(-6*y4^2+6*H*y4)/(H^3-8*y4^3+12*H*y4^2-6*H^2*y4);
    d=(H^3-4*y4^3+9*H*y4^2-6*H^2*y4)/(H^3-8*y4^3+12*H*y4^2-6*H^2*y4);
    
    Aneu(j,1)=a*y(j)^3+b*y(j)^2+c*y(j)+d;
    
end

% Para columna derecha, régimen Poiseuille

for j=1:N
    
    Aneu(j,cfinal)=(-2/H^3)*y(j)^3+(3/H^2)*y(j)^2;
    
end

cquiebre=(R1/R2)*2*(N-1)+1; % número de columnas hasta el quiebre

 % contorno superior completo

for j=2:cquiebre
    Aneu(fh,j)=1;
end

for i=1:fh
    Aneu(i,cquiebre)=1;
end

for j=cquiebre+1:cfinal-1
    Aneu(1,j)=1;
end

for j=2:cfinal-1
    Aneu(N,j)=0.5;
end

% Empieza la iteración

for k=1:2500
    
    
    for i=fh+1:N-1
        for j=2:cquiebre
            Aneu(i,j)=(Aneu(i+1,j)+Aneu(i-1,j)+Aneu(i,j-1)+Aneu(i,j+1))/4;
        end
    end

    for i=2:N-1
        for j=cquiebre+1:cfinal-1
            Aneu(i,j)=(Aneu(i+1,j)+Aneu(i-1,j)+Aneu(i,j-1)+Aneu(i,j+1))/4;
        end
    end
    
end

%% Se crea la matriz B, vorticidad, 

% Se crea la matriz de eta (vorticidad)

    Bneu=zeros(N,cfinal);

% Columna izquierda

    Bneu(fh:N,1)=-(6*a*y(fh:N)+2*b);

% Columna derecha

    Bneu(1:N,cfinal)=(12/H^3)*y(1:N)-6/H^2;

% Contorno izquierdo superior

    Bneu(fh,2:cquiebre)=3*(Aneu(fh,2:cquiebre)-Aneu(fh+1,2:cquiebre))/deltax^2-Bneu(fh+1,2:cquiebre)/2;

% Contorno superior derecho

    Bneu(1,cquiebre+1:cfinal-1)=3*(Aneu(1,cquiebre+1:cfinal-1)-Aneu(1+1,cquiebre+1:cfinal-1))/deltax^2-Bneu(1+1,cquiebre+1:cfinal-1)/2;

% Columna del medio

    Bneu(2:fh-1,cquiebre)=3*(Aneu(2:fh-1,cquiebre)-Aneu(2:fh-1,cquiebre+1))/deltax^2-Bneu(2:fh-1,cquiebre+1)/2;

% Se asume una distribución de vorticidad inicial
    
for i=fh+1:N-1 % fila 5 hasta 7
    Bvector=zeros(1,cfinal-1);
    Bvector=linspace(Bneu(i,1),Bneu(i,cfinal),cfinal);
    Bneu(i,2:(cfinal-1))=Bvector(2:(cfinal-1));
end

for i=2:fh
    Bvector=zeros(1,(cfinal-1)-(cquiebre));
    Bvector=linspace(Bneu(i,cquiebre),Bneu(i,cfinal),(cfinal)-(cquiebre)+1);
    Bneu(i,(cquiebre):(cfinal))=Bvector;
end

%% Proceso iterativo

% Empieza el proceso iterativo, se calcula la vorticidad en los contornos
% superiores de forma aproximada, se calculan las líneas de corriente en
% los nodos interiores, luego la vorticidad en los nodos interiores, luego
% se recalculan las vorticidades en los contornos hasta que converga.

iter=0;
error=10;

while error>maxe
%for iter=1:MaxIt
    
    % Contornos vorticidad 

        % Contorno izquierda, superior 

        for j=2:cquiebre % para columnas 2 a 11
            Bneu(fh,j)=(3.5*Aneu(fh,j)-4*Aneu(fh+1,j)+0.5*Aneu(fh+2,j))/deltax^2;
        end

        % Contorno superior

        for j=cquiebre+1:cfinal-1 % columnas 12 a 18
            Bneu(1,j)=(3.5*Aneu(1,j)-4*Aneu(1+1,j)+0.5*Aneu(1+2,j))/deltax^2;
        end

        % Columna del medio

        for i=2:fh-1 % para fila 2 a 3
            c11=cquiebre;
            Bneu(i,c11)=(3.5*Aneu(i,c11)-4*Aneu(i,c11+1)+0.5*Aneu(i,c11+2))/deltax^2;
        end

        % Puntos de discontinuidad

        Bneu(1,cquiebre)=(Bneu(1,cquiebre+1)+Bneu(2,cquiebre))/2;
        Bneu(fh,cquiebre)=(Bneu(fh-1,cquiebre)+Bneu(fh,cquiebre-1))/2;

    % Se calculan las líneas de corriente interiores

        FS=0.1;

        % Parte central izquierda

        for i=fh+1:N-1 % fila 5 a 7
            for j=2:cquiebre % columna 2 a 11
                Aneu(i,j)=(1-FS)*Aneu(i,j)+(FS/4)*(Aneu(i+1,j)+Aneu(i-1,j)+Aneu(i,j-1)+Aneu(i,j+1)+(deltax^2)*Bneu(i,j));
            end
        end

        % Parte central derecha

        for i=2:N-1 % fila 2 a 7
            for j=cquiebre+1:cfinal-1 % columna 12 a 18
                Aneu(i,j)=(1-FS)*Aneu(i,j)+(FS/4)*(Aneu(i+1,j)+Aneu(i-1,j)+Aneu(i,j-1)+Aneu(i,j+1)+(deltax^2)*Bneu(i,j));
            end
        end
    
	% Vorticidad

        Bneu1=Bneu;

        FV=0.1;
        
        % Parte central izquierda

        for i=fh+1:N-1 % fila 5 a 7
            for j=2:cquiebre % columna 2 a 11

                Bder=Bneu(i,j+1);
                Bizq=Bneu(i,j-1);
                Bsup=Bneu(i-1,j);
                Binf=Bneu(i+1,j);

                beta=(Aneu(i-1,j)-Aneu(i+1,j));
                alfa=(Aneu(i,j+1)-Aneu(i,j-1));

                % Esquema centrado
                
                Bneu(i,j)=(1-FV)*Bneu(i,j)+(FV/4)*((1-(Re/4)*beta)*Bder+(1+(Re/4)*beta)*Bizq+(1+(Re/4)*alfa)*Bsup+(1-(Re/4)*alfa)*Binf);
               
            end
        end

        % Parte central derecha

        for i=2:N-1 % fila 2 hasta 7
            for j=cquiebre+1:cfinal-1 % columna 12 a 18

                Bder=Bneu(i,j+1);
                Bizq=Bneu(i,j-1);
                Bsup=Bneu(i-1,j);
                Binf=Bneu(i+1,j);

                beta=(Aneu(i-1,j)-Aneu(i+1,j));
                alfa=(Aneu(i,j+1)-Aneu(i,j-1));
   
                % Esquema centrado
                
                Bneu(i,j)=(1-FV)*Bneu(i,j)+(FV/4)*((1-(Re/4)*beta)*Bder+(1+(Re/4)*beta)*Bizq+(1+(Re/4)*alfa)*Bsup+(1-(Re/4)*alfa)*Binf);

            end
        end
    
    % Error
    
        error=max(max(abs(Bneu1-Bneu)))
        iter=iter+1
        
        %if iter > 10
        %    error = max(max(abs(Bneu1 - Bneu)))
        %    iter
        %    if error < maxe
        %        break;
        %    end
        %end
        
end
    
for s=1:iteraciones
    
    N2=2*N-1;
    cfinal2=(R1/R2+R3/R2)*2*(N2-1)+1;
    Aneu2=zeros(N2,cfinal2);
    Bneu2=zeros(N2,cfinal2);

    for i=1:N
        for j=1:cfinal
            Aneu2(2*i-1,2*j-1)=Aneu(i,j);
            Bneu2(2*i-1,2*j-1)=Bneu(i,j);
        end
    end

    for i=fh:N-1
        Aneu2(2*i,1:(R1/R2)*2*(N2-1))=(Aneu2(2*i-1,1:(R1/R2)*2*(N2-1))+Aneu2(2*i+1,1:(R1/R2)*2*(N2-1)))/2;
        Bneu2(2*i,1:(R1/R2)*2*(N2-1))=(Bneu2(2*i-1,1:(R1/R2)*2*(N2-1))+Bneu2(2*i+1,1:(R1/R2)*2*(N2-1)))/2;
    end
    
    for i=1:N-1
        Aneu2(2*i,(R1/R2)*2*(N2-1):cfinal2)=(Aneu2(2*i-1,(R1/R2)*2*(N2-1):cfinal2)+Aneu2(2*i+1,(R1/R2)*2*(N2-1):cfinal2))/2;
        Bneu2(2*i,(R1/R2)*2*(N2-1):cfinal2)=(Bneu2(2*i-1,(R1/R2)*2*(N2-1):cfinal2)+Bneu2(2*i+1,(R1/R2)*2*(N2-1):cfinal2))/2;
    end
    
    for j=cquiebre:cfinal-1
        Aneu2(:,2*j)=(Aneu2(:,2*j+1)+Aneu2(:,2*j-1))/2;
        Bneu2(:,2*j)=(Bneu2(:,2*j+1)+Bneu2(:,2*j-1))/2;
    end
    
    for j=1:cquiebre
        Aneu2(ceil(2*(N2-1)*(0.5-1/(2*R2))+0.5):N2,2*j)=(Aneu2(ceil(2*(N2-1)*(0.5-1/(2*R2))+0.5):N2,2*j+1)+Aneu2(ceil(2*(N2-1)*(0.5-1/(2*R2))+0.5):N2,2*j-1))/2;
        Bneu2(ceil(2*(N2-1)*(0.5-1/(2*R2))+0.5):N2,2*j)=(Bneu2(ceil(2*(N2-1)*(0.5-1/(2*R2))+0.5):N2,2*j+1)+Bneu2(ceil(2*(N2-1)*(0.5-1/(2*R2))+0.5):N2,2*j-1))/2;
    end
    
    N=N2;
    cfinal=cfinal2;
    y=linspace(H,H/2,N);
    fh=ceil(2*(N-1)*(0.5-1/(2*R2))+0.5);
    H=2;
    deltax=(H/2)/(N-1); 
    cquiebre=(R1/R2)*2*(N-1)+1;

    Aneu=Aneu2;
    Bneu=Bneu2;
    
    iter=0;
    error=10;
    
    while error>maxe
        %for iter=1:MaxIt

        % Contornos vorticidad 
        
            % Contorno izquierda, superior 

            for j=2:cquiebre % para columnas 2 a 11
                Bneu(fh,j)=(3.5*Aneu(fh,j)-4*Aneu(fh+1,j)+0.5*Aneu(fh+2,j))/deltax^2;
            end

            % Contorno superior

            for j=cquiebre+1:cfinal-1 % columnas 12 a 18
                Bneu(1,j)=(3.5*Aneu(1,j)-4*Aneu(1+1,j)+0.5*Aneu(1+2,j))/deltax^2;
            end

            % Columna del medio

            for i=2:fh-1 % para fila 2 a 3
                c11=cquiebre;
                Bneu(i,c11)=(3.5*Aneu(i,c11)-4*Aneu(i,c11+1)+0.5*Aneu(i,c11+2))/deltax^2;
            end

            % Puntos de discontinuidad

            Bneu(1,cquiebre)=(Bneu(1,cquiebre+1)+Bneu(2,cquiebre))/2;
            Bneu(fh,cquiebre)=(Bneu(fh-1,cquiebre)+Bneu(fh,cquiebre-1))/2;

        % Se calculan las líneas de corriente interiores

            FS=0.1;

            % Parte central izquierda

            for i=fh+1:N-1 % fila 5 a 7
                for j=2:cquiebre % columna 2 a 11
                    Aneu(i,j)=(1-FS)*Aneu(i,j)+(FS/4)*(Aneu(i+1,j)+Aneu(i-1,j)+Aneu(i,j-1)+Aneu(i,j+1)+(deltax^2)*Bneu(i,j));
                end
            end

            % Parte central derecha

            for i=2:N-1 % fila 2 a 7
                for j=cquiebre+1:cfinal-1 % columna 12 a 18
                    Aneu(i,j)=(1-FS)*Aneu(i,j)+(FS/4)*(Aneu(i+1,j)+Aneu(i-1,j)+Aneu(i,j-1)+Aneu(i,j+1)+(deltax^2)*Bneu(i,j));
                end
            end

        % Vorticidad

        Bneu1=Bneu;

        FV=0.1;

            % Parte central izquierda

            for i=fh+1:N-1 % fila 5 a 7
                for j=2:cquiebre % columna 2 a 11

                    Bder=Bneu(i,j+1);
                    Bizq=Bneu(i,j-1);
                    Bsup=Bneu(i-1,j);
                    Binf=Bneu(i+1,j);

                    beta=(Aneu(i-1,j)-Aneu(i+1,j));
                    alfa=(Aneu(i,j+1)-Aneu(i,j-1));

               % Esquema centrado

                    Bneu(i,j)=(1-FV)*Bneu(i,j)+(FV/4)*((1-(Re/4)*beta)*Bder+(1+(Re/4)*beta)*Bizq+(1+(Re/4)*alfa)*Bsup+(1-(Re/4)*alfa)*Binf);

                end
            end

            % Parte central derecha

            for i=2:N-1 % fila 2 hasta 7
                for j=cquiebre+1:cfinal-1 % columna 12 a 18

                    Bder=Bneu(i,j+1);
                    Bizq=Bneu(i,j-1);
                    Bsup=Bneu(i-1,j);
                    Binf=Bneu(i+1,j);

                    beta=(Aneu(i-1,j)-Aneu(i+1,j));
                    alfa=(Aneu(i,j+1)-Aneu(i,j-1));

                    % Esquema centrado

                    Bneu(i,j)=(1-FV)*Bneu(i,j)+(FV/4)*((1-(Re/4)*beta)*Bder+(1+(Re/4)*beta)*Bizq+(1+(Re/4)*alfa)*Bsup+(1-(Re/4)*alfa)*Binf);

                end
            end

            % Error

            error=max(max(abs(Bneu1-Bneu)))
            iter=iter+1

                %if iter > 10
                %    error = max(max(abs(Bneu1 - Bneu)))
                %    iter
                %    if error < maxe
                %        break;
                %    end
                %end

     end
end
    
%% Valores de U, V 

% Se crea la matriz con los valores de velocidad V

V=zeros(N,cfinal);

for i=fh+1:N-1 % fila 5 hasta 7
    for j=2:cquiebre % columna 2 hasta 11
        V(i,j)=(Aneu(i,j+1)-Aneu(i,j-1))/(2*deltax);
    end
end

for i=2:N-1 % desde fila 2 hasta 7
    for j=cquiebre+1:cfinal-1 % desde columna 12 hasta 18
        V(i,j)=(Aneu(i,j+1)-Aneu(i,j-1))/(2*deltax);
    end
end

V=-V;

U=zeros(N,cfinal);

for i=fh+1:N-1 % fila 5 hasta 7
    for j=2:cquiebre % columna 2 hasta 11
        U(fh,j)=-(Aneu(fh,j)-Aneu(fh+1,j))/deltax;
        U(i,j)=(Aneu(i+1,j)-Aneu(i-1,j))/(2*deltax);
        U(N,j)=-(Aneu(N-1,j)-Aneu(N,j))/deltax;
    end
end

% columna izquierda
for i=fh:N
    U(i,1)=-(3*a*y(i)^2+2*b*y(i)+c);
end

% columna derecha

for i=1:N
    U(i,cfinal)=(6/H^3)*y(i)^2-(6/H^2)*y(i);
end

for i=2:N-1 % desde fila 2 hasta 7
    for j=cquiebre+1:cfinal-1 % desde columna 12 hasta 18
        U(1,j)=-(Aneu(1,j)-Aneu(1+1,j))/deltax;
        U(i,j)=(Aneu(i+1,j)-Aneu(i-1,j))/(2*deltax);
        U(N,j)=(Aneu(N,j)-Aneu(N-1,j))/deltax;
    end
end

U=-U;

%% Matriz presión

p=zeros(N,cfinal);

for l=1:5000
    
    p1=p;
    
    for i=fh+1:N-1
        for j=2:cquiebre 
            p(i,j)=((p(i,j+1)+p(i,j-1)+p(i+1,j)+p(i-1,j))-0.5*((U(i,j+1)-U(i,j-1))*(V(i-1,j)-V(i+1,j))-(U(i-1,j)-U(i+1,j))*(V(i,j+1)-V(i,j-1))))/4;
        end
    end

    for i=2:N-1 
        for j=cquiebre+1:cfinal-1
            p(i,j)=((p(i,j+1)+p(i,j-1)+p(i+1,j)+p(i-1,j))-0.5*((U(i,j+1)-U(i,j-1))*(V(i-1,j)-V(i+1,j))-(U(i-1,j)-U(i+1,j))*(V(i,j+1)-V(i,j-1))))/4;
        end
    end
    
    for j=2:cfinal-1 % Fila inferior
        %p(N,j)=p(N,j+1)+0.5*U(N,j)*(U(N,j+1)-U(N,j-1))-(1/Re)*(U(N,j+1)-2*U(N,j)+U(N,j-1)+U(N-2,j)-2*U(N-1,j)+U(N,j))/deltax; %ap/ax
        p(N,j)=p(N-1,j)-(1/Re)*(V(N-2,j)-2*V(N-1,j)+V(N,j))/deltax; %ap/ay
    end
    
    for j=2:cquiebre-1 %Fila izquierda superior
        %p(fh,j)=p(fh,j+1)-(1/Re)*(U(fh,j)-2*U(fh+1,j)+U(fh+2,j))/deltax; %ap/ax
        p(fh,j)=p(fh+1,j)+(1/Re)*(V(fh,j)-2*V(fh+1,j)+V(fh+2,j))/deltax; %ap/ay
    end 
    
    for j=cquiebre+1:cfinal-1 % Fila superior derecha
        %p(1,j)=p(1,j+1)-(1/Re)*(U(1,j)-2*U(2,j)+U(3,j))/deltax; %ap/ax
        p(1,j)=p(1+1,j)+(1/Re)*(V(1,j)-2*V(1+1,j)+V(1+2,j))/deltax; %ap/ay
    end
    
    for i=2:fh-1 %columna quiebre
        %p(i,cquiebre)=p(i-1,cquiebre)-(1/Re)*(V(i,cquiebre)-2*V(i,cquiebre+1)+V(i,cquiebre+2))/deltax; %ap/ay
        p(i,cquiebre)=p(i,cquiebre+1)-(1/Re)*(U(i,cquiebre)-2*U(i,cquiebre+1)+U(i,cquiebre+2))/deltax; %ap/ax
    end
    
    for i=2:N-1 % columna derecha
        p(i,cfinal)=p(i,cfinal-1)+(1/Re)*(U(i-1,cfinal)-2*U(i,cfinal)+U(i+1,cfinal))/deltax; %ap/ax
        %p(i,cfinal)=p(i-1,cfinal)-(1/Re)*(V(i,cfinal-2)-2*V(i,cfinal-1)+V(i,cfinal))/deltax; %ap/ay
        %p(i,cfinal)=0;
    end
    
    for i=fh+1:N-1 % columna izquierda
        p(i,1)=p(i,2)+U(i,1)*(U(i,2)-U(i,1))-(1/Re)*(U(i,1)-2*U(i,2)+U(i,3)+U(i-1,1)-2*U(i,1)+U(i+1,1))/deltax; %ap/ax
        %p(i,1)=p(i-1,1)+U(i,1)*(V(i,2)-V(i,1))-(1/Re)*(V(i,1)-2*V(i,2)+V(i,3))/deltax; %ap/ay
        %p(i,1)=2;
    end
    
    p(fh,1)=1;
    %p(fh,1)=(p(fh,2)+p(fh+1,1))/2;
    p(N,1)=(p(N-1,1)+p(N,2))/2;
    p(fh,cquiebre)=(p(fh,cquiebre-1)+p(fh-1,cquiebre))/2;
    p(1,cquiebre)=(p(2,cquiebre)+p(1,cquiebre+1))/2;
    p(1,cfinal)=(p(1,cfinal-1)+p(2,cfinal))/2;
    p(N,cfinal)=(p(N-1,cfinal)+p(N,cfinal-1))/2;
    
    % Error
        
    if l > 10
        error2 = max(max(abs(p1 - p)))
        l
        if error < maxe
            break;
        end
    end
    
end

figure(2);

[x,y]=meshgrid(0:deltax:(R1+R3)*(h),H/2:deltax:H);
y=flip(y,1);

contour(x,y,p,500)
set(gca,'DataAspectRatio',[1 1 1])

title (['Solución función corriente-vorticidad, presiones, N = ',num2str(N), ', L1/h = ',num2str(R1),', L2/h = ',num2str(R3), ', H/h = ', num2str(R2), ', Re = ', num2str(Re)])

xlabel (' x/h ')
ylabel (' y/h ')

plot(x(N,:),p(N,:))

%% Se grafica

Uneu=flipud(U);
Uuneu=zeros(2*N-1,cfinal);
for i=1:N
    Uuneu(i,:)=U(i,:);
    Uuneu(N+(i-1),:)=Uneu(i,:);
end

Vneu=-flipud(V);
Vvneu=zeros(2*N-1,cfinal);
for i=1:N
    Vvneu(i,:)=V(i,:);
    Vvneu(N+(i-1),:)=Vneu(i,:);
end

%pneu=flipud(p);
%ppneu=zeros(2*N-1,cfinal);
%for i=1:N
%    ppneu(i,:)=p(i,:);
%    ppneu(N+(i-1),:)=pneu(N,:)+pneu(i,:);
%end

[x,y]=meshgrid(0:deltax:(R1+R3)*(h),0:deltax:H);
y=flip(y,1);
figure
%quiver(x,y,Uuneu,Vvneu)

% Líneas de corriente 

figure(1);

plot(0:deltax:(R1+R3)*h,(H/2)*ones(cfinal))
hold on

starty=H/2:deltax:y(fh);
startx=zeros(size(starty));
hlines=streamslice(x,y,Uuneu,Vvneu);
%streamline(x,y,U,V,starty,startx)
hold on

% ploteo contorno

hold on
plot(0:deltax:R1*h,y(fh)*ones(1,R1/R2*2*(N-1)+1),'r')
hold on
plot(0:deltax:R1*h,(H-y(fh))*ones(1,R1/R2*2*(N-1)+1),'r')
hold on
X1=x(1,cquiebre);
Y1(1)=y(fh);
Y1(2)=y(1);
Y2(1)=H-y(fh);
Y2(2)=0;
plot([X1 X1],Y1,'r')
plot([X1 X1],Y2,'r')
hold on
plot(R1*h:deltax:(R1+R3)*h,y(1,1)*ones(1,R3/R2*2*(N-1)+1),'r')
hold on
plot(R1*h:deltax:(R1+R3)*h,zeros(1,R3/R2*2*(N-1)+1),'r')

title (['Solución función corriente-vorticidad, N = ',num2str(N), ', L1/h = ',num2str(R1),', L2/h = ',num2str(R3), ', H/h = ', num2str(R2), ', Re = ', num2str(Re)])

xlabel (' x/h ')
ylabel (' y/h ')

%axis tight

hold off



%figure(3);

%for j=1:10
%    plot(Uuneu(:,j*fix(cfinal/10)),y(:,1))
%    hold on
%end
%hold off
    
%figure (4);
%for j=1:cfinal
%    plot(ppneu(:,j),y(:,1))
%    hold on
%end
%hold off

% Líneas de corriente 
