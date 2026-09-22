%% Función vorticidad - líneas de corriente, expansión abrupta de un flujo, caso axialsimetrico

%%

% Se definen valores 

H=2;
N=113; % número de nodos
dr=(H/2)/(N-1); 
Re=100;
MaxIt=350000;
maxe = 5e-7;

% Relación H/h

R2=2; 
h=H/R2;

% Relación L1/h

R1=1;

% Relación L2/h

R3=8;

% Se crea la matriz de líneas de corriente

cfinal=(R1/R2+R3/R2)*2*(N-1)+1; % número de columnas 

%% Se crea la matriz A, líneas de corriente, con las condiciones de contorno y se resuelve el caso potencial

Aneu=zeros(N,cfinal);

% Se agregan las condiciones de contorno 

% Columna izquierda, régimen Poiseuille

r=linspace(H,H/2,N);

fh=ceil(2*(N-1)*(0.5-1/(2*R2))+0.5); % fila del valor de h 

for j=fh:N % para la fila 4 hasta 8

    r4=r(fh);

    a=2/(H^3-8*r4^3+12*H*r4^2-6*H^2*r4);
    b=(-3*H)/(H^3-8*r4^3+12*H*r4^2-6*H^2*r4);
    c=(-6*r4^2+6*H*r4)/(H^3-8*r4^3+12*H*r4^2-6*H^2*r4);
    d=(H^3-4*r4^3+9*H*r4^2-6*H^2*r4)/(H^3-8*r4^3+12*H*r4^2-6*H^2*r4);
    
    Aneu(j,1)=a*r(j)^3+b*r(j)^2+c*r(j)+d;
    
end

% Para columna derecha, régimen Poiseuille

for j=1:N
    
    Aneu(j,cfinal)=(-2/H^3)*r(j)^3+(3/H^2)*r(j)^2;
    
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
            Aneu(i,j)=(Aneu(i+1,j)+Aneu(i-1,j)+Aneu(i,j-1)+Aneu(i,j+1)-dr/(2*r(i))*(Aneu(i-1,j)-Aneu(i+1,j)))/4;
        end
    end

    for i=2:N-1
        for j=cquiebre+1:cfinal-1
            Aneu(i,j)=(Aneu(i+1,j)+Aneu(i-1,j)+Aneu(i,j-1)+Aneu(i,j+1)-dr/(2*r(i))*(Aneu(i-1,j)-Aneu(i+1,j)))/4;
        end
    end
    
end

%% Se crea la matriz B, vorticidad, 

% Se crea la matriz de eta (vorticidad)

    Bneu=zeros(N,cfinal);

% Columna izquierda

    Bneu(fh:N,1)=-(6*a*r(fh:N)+2*b);

% Columna derecha

    Bneu(1:N,cfinal)=(12/H^3)*r(1:N)-6/H^2;

% Contorno izquierdo superior

    Bneu(fh,2:cquiebre)=3*(Aneu(fh,2:cquiebre)-Aneu(fh+1,2:cquiebre))/dr^2-Bneu(fh+1,2:cquiebre)/2;

% Contorno superior derecho

    Bneu(1,cquiebre+1:cfinal-1)=3*(Aneu(1,cquiebre+1:cfinal-1)-Aneu(1+1,cquiebre+1:cfinal-1))/dr^2-Bneu(1+1,cquiebre+1:cfinal-1)/2;

% Columna del medio

    Bneu(2:fh-1,cquiebre)=3*(Aneu(2:fh-1,cquiebre)-Aneu(2:fh-1,cquiebre+1))/dr^2-Bneu(2:fh-1,cquiebre+1)/2;

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


for iter=1:MaxIt
    
    % Contornos vorticidad 

        % Contorno izquierda, superior 

        for j=2:cquiebre % para columnas 2 a 11
            %Bneu(fh,j)=2*(Aneu(fh,j)-Aneu(fh+1,j))/(r(j)*dr^2)-(2/r(fh))*(Aneu(fh,j)-Aneu(fh+1,j))/dr^2+(1/r(fh)^2))*(Aneu(fh,j)-Aneu(fh+1,j))/dr;
            %Bneu(fh,j)=(-(Aneu(fh,j)-2*Aneu(fh+1,j)+Aneu(fh+2,j))/dr^2+(1/r(i,j))*(Aneu(fh,j)-Aneu(fh+1,j)))/r(i,j);
            Bneu(fh,j)=2*(Aneu(fh,j)-Aneu(fh+1,j))/(r(fh)*dr^2);
        end

        % Contorno derecho, superior

        for j=cquiebre+1:cfinal-1 % columnas 12 a 18
            Bneu(1,j)=2*(Aneu(1,j)-Aneu(1+1,j))/(r(1)*dr^2);
            %Bneu(fh,j)=2*(Aneu(fh,j)-Aneu(fh+1,j))/(r(fh)*dr^2)-(2/r(j))*(Aneu(fh,j)-Aneu(fh+1,j))/dr^2+(1/r(fh)^2))*(Aneu(fh,j)-Aneu(fh+1,j))/dr;
        end

        % Columna de quiebre

        for i=2:fh-1 % para fila 2 a 3
            Bneu(i,cquiebre)=2*(Aneu(i,cquiebre)-Aneu(i,cquiebre+1))/(r(i)*dr^2);
        end

        % Puntos de discontinuidad

        Bneu(1,cquiebre)=(Bneu(1,cquiebre+1)+Bneu(2,cquiebre))/2;
        Bneu(fh,cquiebre)=(Bneu(fh-1,cquiebre)+Bneu(fh,cquiebre-1))/2;

    % Se calculan las líneas de corriente interiores

        FS=0.05;

        % Parte central izquierda

        for i=fh+1:N-1 % fila 5 a 7
            for j=2:cquiebre % columna 2 a 11
                Aneu(i,j)=(1-FS)*Aneu(i,j)+0.25*FS*(Aneu(i+1,j)+Aneu(i-1,j)+Aneu(i,j-1)+Aneu(i,j+1)-dr/(2*r(i))*(Aneu(i-1,j)-Aneu(i+1,j))+Bneu(i,j)*r(i)*dr^2);
            end
        end

        % Parte central derecha

        for i=2:N-1 % fila 2 a 7
            for j=cquiebre+1:cfinal-1 % columna 12 a 18
                Aneu(i,j)=(1-FS)*Aneu(i,j)+0.25*FS*(Aneu(i+1,j)+Aneu(i-1,j)+Aneu(i,j-1)+Aneu(i,j+1)-dr/(2*r(i))*(Aneu(i-1,j)-Aneu(i+1,j))+Bneu(i,j)*r(i)*dr^2);
            end
        end
    
	% Vorticidad

        Bneu1=Bneu;

        FV=0.05;
        
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
                
                Bneu(i,j)=(1-FV)*Bneu(i,j)+FV*0.25*((Bder+Bizq+Bsup+Binf)+(dr/r(i))*(Bsup-Binf)-(Re/(4*r(i)))*beta*(Bder-Bizq)+(Re/(4*r(i)))*alfa*(Bsup-Binf));
               
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
                
                Bneu(i,j)=(1-FV)*Bneu(i,j)+FV*0.25*((Bder+Bizq+Bsup+Binf)+(dr/r(i))*(Bsup-Binf)-(Re/(4*r(i)))*beta*(Bder-Bizq)+(Re/(4*r(i)))*alfa*(Bsup-Binf));
               
            end
        end
    
    % Error
        
        if iter > 10
            error = max(max(abs(Bneu1 - Bneu)))
            iter
            if error < maxe
                break;
            end
        end
        
end
    
%% Valores de U, V 

% Se crea la matriz con los valores de velocidad V

V=zeros(N,cfinal);

for i=fh+1:N-1 % fila 5 hasta 7
    for j=2:cquiebre % columna 2 hasta 11
        V(i,j)=(Aneu(i,j+1)-Aneu(i,j-1))/(2*dr);
    end
end

for i=2:N-1 % desde fila 2 hasta 7
    for j=cquiebre+1:cfinal-1 % desde columna 12 hasta 18
        V(i,j)=(Aneu(i,j+1)-Aneu(i,j-1))/(2*dr);
    end
end

V=-V;

U=zeros(N,cfinal);

for i=fh+1:N-1 % fila 5 hasta 7
    for j=2:cquiebre % columna 2 hasta 11
        U(fh,j)=-(Aneu(fh,j)-Aneu(fh+1,j))/dr;
        U(i,j)=(Aneu(i+1,j)-Aneu(i-1,j))/(2*dr);
        U(N,j)=-(Aneu(N-1,j)-Aneu(N,j))/dr;
    end
end

% columna izquierda
for i=fh:N
    U(i,1)=-(3*a*r(i)^2+2*b*r(i)+c);
end

% columna derecha

for i=1:N
    U(i,cfinal)=(6/H^3)*r(i)^2-(6/H^2)*r(i);
end

for i=2:N-1 % desde fila 2 hasta 7
    for j=cquiebre+1:cfinal-1 % desde columna 12 hasta 18
        U(1,j)=-(Aneu(1,j)-Aneu(1+1,j))/dr;
        U(i,j)=(Aneu(i+1,j)-Aneu(i-1,j))/(2*dr);
        U(N,j)=(Aneu(N,j)-Aneu(N-1,j))/dr;
    end
end

U=-U;

%% Matriz presión

p=zeros(N,cfinal);

for l=1:10000
    
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
        p(N,j)=p(N-1,j)-(1/Re)*(V(N-2,j)-2*V(N-1,j)+V(N,j))/dr; %ap/ay
    end
    
    for j=2:cquiebre-1 %Fila izquierda superior
        %p(fh,j)=p(fh,j+1)-(1/Re)*(U(fh,j)-2*U(fh+1,j)+U(fh+2,j))/deltax; %ap/ax
        p(fh,j)=p(fh+1,j)+(1/Re)*(V(fh,j)-2*V(fh+1,j)+V(fh+2,j))/dr; %ap/ay
    end 
    
    for j=cquiebre+1:cfinal-1 % Fila superior derecha
        %p(1,j)=p(1,j+1)-(1/Re)*(U(1,j)-2*U(2,j)+U(3,j))/deltax; %ap/ax
        p(1,j)=p(1+1,j)+(1/Re)*(V(1,j)-2*V(1+1,j)+V(1+2,j))/dr; %ap/ay
    end
    
    for i=2:fh-1 %columna quiebre
        %p(i,cquiebre)=p(i-1,cquiebre)-(1/Re)*(V(i,cquiebre)-2*V(i,cquiebre+1)+V(i,cquiebre+2))/deltax; %ap/ay
        p(i,cquiebre)=p(i,cquiebre+1)-(1/Re)*(U(i,cquiebre)-2*U(i,cquiebre+1)+U(i,cquiebre+2))/dr; %ap/ax
    end
    
    for i=2:N-1 % columna derecha
        p(i,cfinal)=p(i,cfinal-1)+(1/Re)*(U(i-1,cfinal)-2*U(i,cfinal)+U(i+1,cfinal))/dr; %ap/ax
        %p(i,cfinal)=p(i-1,cfinal)-(1/Re)*(V(i,cfinal-2)-2*V(i,cfinal-1)+V(i,cfinal))/deltax; %ap/ay
        %p(i,cfinal)=0;
    end
    
    for i=fh+1:N-1 % columna izquierda
        p(i,1)=p(i,2)+U(i,1)*(U(i,2)-U(i,1))-(1/Re)*(U(i,1)-2*U(i,2)+U(i,3)+U(i-1,1)-2*U(i,1)+U(i+1,1))/dr; %ap/ax
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

[x,r]=meshgrid(0:dr:(R1+R3)*(h),H/2:dr:H);
r=flip(r,1);

contour(x,r,p,500)
set(gca,'DataAspectRatio',[1 1 1])

title (['Solución función corriente-vorticidad, presiones, N = ',num2str(N), ', L1/h = ',num2str(R1),', L2/h = ',num2str(R3), ', H/h = ', num2str(R2), ', Re = ', num2str(Re)])

xlabel (' x/h ')
ylabel (' y/h ')

%plot(x(N,:),p(N,:))

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

[x,r]=meshgrid(0:dr:(R1+R3)*(h),0:dr:H);
r=flip(r,1);
figure
%quiver(x,y,Uuneu,Vvneu)

% Líneas de corriente 

figure(1);

plot(0:dr:(R1+R3)*h,(H/2)*ones(cfinal))
hold on

starty=H/2:dr:r(fh);
startx=zeros(size(starty));
hlines=streamslice(x,r,Uuneu,Vvneu);
%streamline(x,y,U,V,starty,startx)
hold on

% ploteo contorno

hold on
plot(0:dr:R1*h,r(fh)*ones(1,R1/R2*2*(N-1)+1),'r')
hold on
plot(0:dr:R1*h,(H-r(fh))*ones(1,R1/R2*2*(N-1)+1),'r')
hold on
X1=x(1,cquiebre);
Y1(1)=r(fh);
Y1(2)=r(1);
Y2(1)=H-r(fh);
Y2(2)=0;
plot([X1 X1],Y1,'r')
plot([X1 X1],Y2,'r')
hold on
plot(R1*h:dr:(R1+R3)*h,r(1,1)*ones(1,R3/R2*2*(N-1)+1),'r')
hold on
plot(R1*h:dr:(R1+R3)*h,zeros(1,R3/R2*2*(N-1)+1),'r')

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
