  
%truss analysis
clearvars;clc;
%input
COORD =[0 0; 10 0; 20 0; 0 12; 10 12; 20 12]; %node coordinate. node 1 is in center point
CON = [1 2; 2 3; 4 5; 5 6; 1 4; 2 4; 1 5; 2 5; 3 5; 2 6; 3 6]; %element conncting node
EQ = [10 11; 1 2; 3 12; 4 5; 6 7; 8 9];%each node reaction force. Degree of freedom variable p1,p2---p12
NR = 3;%no of support reaction force
NE = size(CON,1);%no of elements
NN = size(COORD,1); %no of nodes
EA = [2 3 3 3 2 2 3 2 3 4 2]'*10^3; %E * A value are different for each elements
Pf = [0 0 0 0 0 0 -10 0 0]'; %external load from DOF p1 to p9. p10,p11,p12 are fixed point
Ur = [0 0 0]';% P10,P11, P12 are the support reaction

%calculation
%Structural information
NOS = NE+NR-2*NN;%no of static indetermanacy
NOK = 2*NN-NR;%no of kinematic indetermanacy

%length of the elements
L = zeros(NE,1);
for k =1:NE
    i= CON(k,1);%first starting point of local node
    j=CON(k,2);%ending point of local node
    dx=COORD(j,1)-COORD(i,1);
    dy=COORD(j,2)-COORD(i,2);
    L(k) =sqrt(dx^2+dy^2);
end

%degree of freedon variable(p1,p2,-----p12)-ID array
ID= zeros(NE,4);
for k = 1:NE
     i= CON(k,1);%first starting point of local node
     j=CON(k,2);%ending point of local node
     ID(k,1:2)=EQ(i,1:2);
     ID(k,3:4)=EQ(j,1:2);
     
end


%stiffness matrix
NDOF= 2*NN;
K=zeros (NDOF, NDOF);
for k =1:NE
    i=CON(k,1);%first starting point of local node
    j=CON(k,2);%ending point of local node
    dx=COORD(j,1) - COORD(i,1);
    dy=COORD(j,2) - COORD(i,2);
    a= [-dx/L(k) -dy/L(k) dx/L(k) dy/L(k)];%cos(theta)=dx/L(k)),sin(theta)=dy/L(k) 
    ES = a' .*EA(k)/L(k)*a;
    %assembly of global stiffness matrix
    for m = 1:4
        for n =1:4
            mi =ID(k,m);
            ni = ID(k,n);
            K(mi,ni)= K(mi,ni) + ES(m,n);
        end
    end
       
end
fprintf('Global stiffness matrix K')
K
Kff(1:NOK,1:NOK) =  K(1:NOK,1:NOK);
Kfr(1:NOK,1:NDOF-NOK) = K(1:NOK, NOK+1:NDOF);
Krf=Kfr';
Krr(1:NDOF-NOK,1:NDOF-NOK) =K(NOK+1:NDOF,NOK+1:NDOF);


%deformation 
Uf = Kff\Pf;%guess elemination
fprintf('Deflection in each element')
U =[Uf;Ur]%deflection
scale = 10;

%internal force
N = zeros(NE,1);
for  k =1:NE
     i=CON(k,1);%first starting point of local node
     j=CON(k,2);%ending point of local node
     dx=COORD(j,1) - COORD(i,1);
     dy=COORD(j,2) - COORD(i,2);
     a= [-dx/L(k) -dy/L(k) dx/L(k) dy/L(k)];%cos(theta)=dx/L(k)),sin(theta)=dy/L(k)
     u =zeros(4,1);
     for m=1:4
         u(m) = U(ID(k,m));
     end
     N(k) = EA(k)/L(k).*a*u;
end
%support reaction
R=Krf*Uf + Krr*Ur
%----------------------------------------------------------------------
%Plot structure
f1=figure();
NCOORD = zeros(size(COORD));%deformed co-ordinate generation through zero matrix
scale = 10;
for n =1:NN
    NCOORD(n,1) = COORD(n,1) +scale*U(EQ(n,1));
    NCOORD(n,2) = COORD(n,2) +scale*U(EQ(n,2));
end

for k =1:NE
     i=CON(k,1);%first starting point of local node
     j=CON(k,2);%ending point of local node
     x=[COORD(i,1) COORD(j,1)];
     y=[COORD(i,2) COORD(j,2)];
     xlim([-1 21]);
     ylim([-1 13]);
     plot(x,y,'k-');
     hold on
     ux=[NCOORD(i,1) NCOORD(j,1)];
     uy=[NCOORD(i,2) NCOORD(j,2)];
     xlim([-1 21]);
     ylim([-1 13]);
     plot(ux,uy,'r--');
     hold on
end




