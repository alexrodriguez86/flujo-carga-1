# flujo-carga
 flujo de carga NEWTON RAPHSON
 %FLUJO DE CARGA MODIFICADO POR PULL REQUEST
 %CARLOS SALDARRIAGA UNIVERSIDAD TECNOLOGICA DE PEREIRA
 %FECHA: 7 FEB 2019

 clc        %limpia pantalla
 clear all  %borra variables
 close all  %cierra ventanas

 %carga de datos

 %nombre archivo
 [LIBRO DIR] = uigetfile('*.xlsx','archivo excel');

 %letrero espera

 wai = waitbar(0,'por favor espere');

 %tipo  1>>slack;  2>>PV;  3>>PQ

 %Datos nodos

 %Nodos = [bus_i tipo Pg Qg Pd Qd Bsh Vm Va baseKV Vmax Vmin]
 Nodos = xlsread(LIBRO,'Nodos');

 waitbar(0.5,wai)

 %datos entre lineas
 %Lineas = [Ni Nf R X Bsh/2]
 Lineas = xlsread(LIBRO,'Lineas');

 waitbar(1,wai)
 close(wai)

 NR=size(Lineas,1);%numero de lineas
 NB=size(Nodos,1) ;%numero de barras o nodos

 %impedancia serie de lineas
 Zlinea=Lineas(:,3)+(sqrt(-1)*Lineas(:,4));

 %valores iniciales
 Vm_o=Nodos(:,8);%magnitud de las tensiones
 Va_o=Nodos(:,9);%angulo de las tensiones

 %tolerancia
 Tol=0.001;

 %construccion de matriz de admitancia nodal (ybus)

 %elementos de la diagonal
 for k = 1:NR

     Ybij=-1/Zlinea(k);%admitancia serie nodo inicio y fin
     i=Lineas(k,1);%nodo de inicio
     j=Lineas(k,2);%nodo fin
     Ybus(i,j)=Ybij;%posicion ij en la YBUS
     Ybus(j,i)=Ybij;%se iguala por ser una matriz simetrica

 end

 %elementos de la diagonal.

 for i = 1:NB%nuemero de barras o nodos
     for j = 1:NR%lineas de transmision

         if i==Lineas(j,1)
             Ybus(i,i)=Ybus(i,i)+(1/Zlinea(j))+(sqrt(-1)*Lineas(j,5));
         end

         if i==Lineas(j,2)
             Ybus(i,i)=Ybus(i,i)+(1/Zlinea(j))+(sqrt(-1)*Lineas(j,5));
         end

     end

     Ybus(i,i) = Ybus(i,i) + Nodos(i,7);

 end

 %figure
 %spy(Ybus)
 %title('estructura de la Ybus')

 %contruye variables magnitud y angulos de los voltajes

 for k = 1:NB
     syms (['V' num2str(k)]);
     Vm(k,1) = eval(['V',num2str(k)]);

     syms (['Teta' num2str(k)]);
     Va(k,1) = eval(['Teta',num2str(k)]);
 end

 i=sqrt(-1);

 %contruye ecuaciones de balance nodal

 for k1=1:NB
     for k2=1:NB
         T(k1,k2)=Va(k2)-Va(k1);

         syms (['V' num2str(k2)]);
         V_(k1,k2) = eval(['V',num2str(k2)]);
     end
 end
 G_=real(Ybus);
 B_=imag(Ybus);

 Dp=(Nodos(:,3)-Nodos(:,5))-((diag((G_.*V_)*cos(T)).*Vm)+(diag((B_.*V_)*sin(T)).*Vm));
 Dq=(Nodos(:,4)-Nodos(:,6))-((diag((G_.*V_)*sin(T)).*Vm)-(diag((B_.*V_)*cos(T)).*Vm));

 %seleccion para los delta Dp (PV,PQ) y para los Dq
 for h=1:NB
     if (Nodos(h,2)==1)%slack
         Dp(h)=[inf];
         Dq(h)=[inf];

         Vm(h)=inf;
         Va(h)=inf;

     end

     if (Nodos(h,2)==2)%PV
         Dq(h)=[inf];

         Vm(h)=inf;

     end

 end

 Dp(find(inf==Dp))=[];
 Dq(find(inf==Dq))=[];

 Vm(find(inf==Vm))=[];
 Va(find(inf==Va))=[];

 X=[Va;Vm];%variables de estado desconocidas
 Y=[Dp;Dq];%sistema de ecuaciones no lineales

 %%calculo del jacobiano [jacobian (X,Y)]

 for k1=1:length(Y)
     for k2=1:length(X)
         J(k1,k2)=diff(Y(k1),X(k2));
     end
 end

 %metodo iterativo de newton raphson

 for h = 1:NB
     evalc(['V' num2str(h) ' = Vm_o(h)']);
     evalc(['Teta' num2str(h) ' = Va_o(h)']);

 end

 %calcula delta P u delta Q

 Delta=eval(Y);
 ITE=1;
 while (ITE<=100)&(sum((abs(Delta)>=Tol))>0)

     J_o=eval(-1*J); %evalua el Jacobiano
     Delt_T_V=inv(J_o)*Delta;%calcula delta teta y delta V

     pv=length(find(2==Nodos(:,2)));%numero de nodos pv
     pq=length(find(3==Nodos(:,2)));%numero de nodos pq

     %actualiza delta teta
     Va_o((find(2<=Nodos(:,2))))=Va_o((find(2<=Nodos(:,2))))+Delt_T_V(1:pv+pq);

     %actualiza delta V
     Vm_o((find(3==Nodos(:,2))))=Vm_o((find(3==Nodos(:,2))))+Delt_T_V(pv+pq+1:length(Delt_T_V));

     for h = 1:NB
         evalc(['V' num2str(h) ' = Vm_o(h)']);
         evalc(['Teta' num2str(h) ' = Va_o(h)']);

     end

     Delta=eval(Y);

     ITE=ITE+1;

 end

 %claculo de flujo por las lineas

 En=(Vm_o.').*(cos((Va_o.'))+i*sin((Va_o.')));%tension nodal fasorial
 A=zeros(NR,NB);
 for k1=1:NR

     A(k1,Lineas(k1,1))=1;
     A(k1,Lineas(k1,2))=-1;

     A1(k1,:)=A(k1,:).*En;
     B1(k1,1)=i*(Lineas(k1,5).*En(Lineas(k1,1)));
     B2(k1,1)=i*(Lineas(k1,5).*En(Lineas(k1,2)));
     Em1(k1,1)=En(Lineas(k1,1));
     Em2(k1,1)=En(Lineas(k1,2));

 end

 Iij=((1./Zlinea).*sum( A1,2))+B1;
 Iji=((1./Zlinea).*sum(-A1,2))+B2;

 Sij_o=Em1.*conj(Iij);
 Sji_o=Em2.*conj(Iji);

 Pij=real(Sij_o);
 Qij=imag(Sij_o);

 Pji=real(Sji_o);
 Qji=imag(Sji_o);

 Si=(En.').*conj(Ybus*(En.'))%potencia aparente nodal
 Pi=real(Si);%potencia activa
 Qi=imag(Si);%potencia reactiva



 %muestra de resultados

 disp('________________________________')
 fprintf('Nï¿½mero de iteraciones: %d\n',ITE)
 disp('________________________________')
 disp('  ')
 disp('________________________________')
 disp('magnitud de las Tensiones nodales pu')
 disp('________________________________')
 for k=1:length(Vm_o)
     fprintf('  V_%d = %.3f',k,Vm_o(k))
     if mod(k,2)==0 % si k es par
         fprintf('\n')

     else
         fprintf('      ;')

     end

 end

 disp('________________________________')
 disp('   ')
 disp('________________________________')
 disp('Angulo de las tensiones nodales (ï¿½)')
 for k=1:length(Va_o)
     if Va_o(k)<0 %angulo negativo
         fprintf('    V_%d = %.2f',k,Va_o(k)*180/pi)
     else
         fprintf('    V_%d = %.3f',k,Va_o(k)*180/pi)

     end

     if mod(k,2)==0 %si k es par
         fprintf('\n')

     else
         fprintf('    ;')

     end

 end

 disp('________________________________')
 disp('   ')
 disp('________________________________')
 disp('Potencias nodales (PU)')
 disp('________________________________')

 for k=1:length(Pi)
     if Pi(k)<0 %Pi negativo
         fprintf('    P_%d = %.2f   ;',k,Pi(k))
     else
         fprintf('    P_%d = %.3f   ;',k,Pi(k))

     end

     if Qi(k)<0 %Qi k es negativo
         fprintf('    Q_%d = %.2f\n',k,Qi(k))
     else
         fprintf('    Q_%d = %.3f\n',k,Qi(k))



     end

 end

 disp('________________________________')
 disp('   ')
 disp('________________________________')
 disp('Potencias por las lineas (PU)')
 disp('________________________________')

 for k=1:length(Pij)
     k1=Lineas(k,1:2);
     if Pij(k)<0 %Pij negativo
         fprintf('    P_%d_%d = %.2f ;',k1(1),k1(2),Pij(k))
     else
         fprintf('    P_%d_%d = %.3f ;',k1(1),k1(2),Pij(k))

     end

     if Qij(k)<0 %Qi k es negativo
         fprintf('    Q_%d_%d = %.2f \n',k1(1),k1(2),Qij(k))
     else
         fprintf('    Q_%d_%d = %.3f \n',k1(1),k1(2),Qij(k))


     end

 end

 disp('__________________________________________')


 for k=1:length(Pij)
     k1=Lineas(k,1:2);
     if Pji(k)<0 %Pi negativo
         fprintf('    P_%d_%d = %.2f ;',k1(2),k1(1),Pji(k))
     else
         fprintf('    P_%d_%d = %.3f ;',k1(2),k1(1),Pji(k))

     end

     if Qji(k)<0 %Qi k es negativo
         fprintf('    Q%d_%d = %.2f \n',k1(2),k1(1),Qji(k))
     else
         fprintf('    Q%d_%d = %.3f \n',k1(2),k1(1),Qji(k))


     end

 end

 disp('__________________________________________')

 %graficas

 %magnitud tensiones

 subplot(1,2,1)
 bar(Vm_o)
 grid on
 title('tension nodal')
 xlabel('nodo')
 ylabel('pu')

 %potencia neta inyectada

 subplot(1,2,2)
 bar([Pi,Qi])
 grid on
 legend('P_i','Q_i')
 title('potencia neta inyectada')
 xlabel('nodo')
 ylabel('pu')
