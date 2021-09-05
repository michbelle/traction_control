%% UKF
% static values

%Q=E(ww')
Q=[...
    1 0 0 0 0 0 0 0;...
    0 2 0 0 0 0 0 0;...
    0 0 3 0 0 0 0 0;...
    0 0 0 4 0 0 0 0;...
    0 0 0 0 5 0 0 0;...
    0 0 0 0 0 6 0 0;...
    0 0 0 0 0 0 7 0;...
    0 0 0 0 0 0 0 8];

R=[...
    1 0 0 0 0 0; ...
    0 2 0 0 0 0; ...
    0 0 3 0 0 0; ...
    0 0 0 4 0 0; ...
    0 0 0 0 5 0; ...
    0 0 0 0 0 6];

%x_0=         [vx; vy; d_phi; w_fl; w_fr; w_rl; w_rr; mu_max];
x_0=          [  0;  0;     0;    0;    0;    0;    0;   1];
predicted_x_0=[  0;  0;     0;    0;    0;    0;    0;   1];
P_x_0=(x_0-predicted_x_0)*(x_0-predicted_x_0)';


%k value for scalling and calcolo sigma points
n=8;
m=6;
L=n+m;
k_ukf=L-3;
%pesi sigma points
w0=k_ukf/(L+k_ukf);
w1=1/(2*L+2*k_ukf);
%2L+1=28 TOCHECK
w=[w0;w1;w1;w1;w1;...
    w1;w1;w1;w1;w1;...
    w1;w1;w1;w1;w1;...
    w1;w1;w1;w1;w1;...
    w1;w1;w1;w1;w1;...
    w1;w1;w1;w1];

%valore scalare distribuzione sigma points

eta_ukf=sqrt(L+k_ukf);

%An efficient way for calculating the square root 
%factorization is given by the Cholesky decomposition.

%TOCHECK
%used integrated matlab fuction

X_w=[...
    [0 0 0 0 0 0 0 0]' ...
    eta_ukf*sqrtm(Q) ...
    -eta_ukf*sqrtm(Q)];%n x 2n+1

G_v=[...
    eta_ukf*sqrtm(R) ...
    -eta_ukf*sqrtm(R)];%mx2m


%% UKF
x_k=x_0;
P_x=P_x_0;

X_x_k=[...
    x_k ...
    x_k+eta_ukf*sqrtm(P_x) ...
    x_k-eta_ukf*sqrtm(P_x)];%n x 2n+1

%propagazione sigma points nel sistema
%i=0..2n 

%insert input
%u_k=[d_fl;d_fr;M_fl;M_fr;M_rl;M_rr];
u_k=[0;1;0;0;0;0;0;0];

%TODO di test dopo riporto

X_final=[];

for i=1:2*n+1
    
    X_temp=(X_x_k(:,i)+u_k)+X_w(:,i);
    X_final=[X_final X_temp]; 
    
end

%extended state sigma points vector to fill until 2*L+1=29 colonne

less=2*L+1-size(X_final,2);
X_k1=[X_final zeros(8,less)];%used in update step

%tochastic properties of the random state variables
%are approximated up to the second order 

x_k1=[];

x_k1=w(1)*X_x_k(:,1);

for i=2:2*n+1
    x_k1_temp=w(i)*X_x_k(:,i);
    x_k1=x_k1+x_k1_temp;
end

P_x_k1=w(1)*(X_x_k(:,1)-x_k1)*(X_x_k(:,1)-x_k1)';

for i=2:2*n+1
    P_x_k1_temp=w(i)*(X_x_k(:,i)-x_k1)*(X_x_k(:,i)-x_k1)';
    P_x_k1=P_x_k1+P_x_k1_temp;
end

%output sigma points 


%mx2L+1

Y_final=[];

for i=1:2*n+1
    
    Y_temp=[0 0 0 0 0 1]';%(X_x_k(:,i))+ u ricava funzione da sigma points
    Y_final=[Y_final Y_temp]; 
    
end

%mx2L+1 = 29
%Y_final=17
%G_v =12

Y_k1=[Y_final Y_final(:,1)+G_v];


%output stima

y_k1=[];

y_k1=w(1)*Y_k1(:,1);

for i=2:2*L+1
    y_k1_temp=w(i)*Y_k1(:,i);
    y_k1=y_k1+y_k1_temp;
end


P_y_k1=w(1)*(Y_k1(:,1)-y_k1)*(Y_k1(:,1)-y_k1)';

for i=2:2*L+1
    P_y_k1_temp=w(i)*(Y_k1(:,i)-y_k1)*(Y_k1(:,i)-y_k1)';
    P_y_k1=P_y_k1+P_y_k1_temp;
end



%Kalman gain


P_xy_k1=w(1)*(X_k1(:,1)-x_k1)*(Y_k1(:,1)-y_k1)';

for i=2:2*L+1
    P_xy_k_temp=w(i)*(X_k1(:,i)-x_k1)*(Y_k1(:,i)-y_k1)';
    P_xy_k1=P_xy_k1+P_xy_k_temp;
end


K_k1=P_xy_k1/(P_y_k1);%/=inv()

%new measure y
y=[0; 0; 0; 0; 0; 0];
x_k1k1=x_k1+K_k1*(y-y_k1);


P_x_k1k1= P_x_k1-K_k1*P_y_k1*K_k1';







