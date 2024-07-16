clc;
clear;
close all;

G=8*(10^10);
Ro=7800;
Dl=0.25;
D0=0.75;
L=4;
Kt=50000000;

x = 0:0.01:L;
Dx=D0-((Dl-D0)/L)*x;
Ipx=(pi*Dx.^4)/32;
c=(G/Ro)^0.5;

%% det
w = [1.6680    4.7461    7.8737   11.0090   14.1469];

F = [];
G = [];

for j = 1:5
    for i = 1:5
        f = @(x) sin(w(i)*x) * sin(w(j)*x);
        g = @(x) (-2/(0.25+0.5*(1-x/4)))*(w(i))*cos(w(i)*x)*sin(w(j)*x);
        F(j,i) = integral(f, 0, 4, 'ArrayValued', true); 
        G(j,i) = integral(g, 0, 4, 'ArrayValued', true);  
    end
end

syms B
for j = 1:5
    for i = 1:5
        W(j, i) = (B - w(i)^2) * F(j,i) + G(j,i);
    end
end

A = det(W);
s = solve(A == 0, B);
s = double(s);
s= sort(s);

disp('%%%%%%%%%%%%%%%%%%%% "w^2/c^2 = B" B_ha:');
disp(' ');

disp(s);
wn = []; 
for i = 1:5
    wn(1,i) = sqrt(s(i,1))*sqrt(8*10^10/7800)/4; 
end
fprintf('natural frequencies of the system are:\n')
disp(wn);

C = zeros(5, 5);

for i = 1:5
    c1 = subs(W, s(i, 1));
    syms C1 C2 C3 C4
    eqn1 = c1(1, :) * [C1; C2; C3; C4; 1] == 0;
    eqn2 = c1(2, :) * [C1; C2; C3; C4; 1] == 0;
    eqn3 = c1(3, :) * [C1; C2; C3; C4; 1] == 0;
    eqn4 = c1(4, :) * [C1; C2; C3; C4; 1] == 0;
    
    solutions = solve(eqn1, eqn2, eqn3, eqn4);
    
    C(:, i) = double([solutions.C1; solutions.C2; solutions.C3; solutions.C4 ; 1]);
    C(:,i) = C(:,i)/(max(abs(C(:,i))));
end

% for i=1:5
%     C(:,i) = C(:,i)/(max(abs(C(:,i))));
% end

si = [];
x = 0:0.01:4;
for i = 1:5
    si(i,:) = [sin(wn(1,i)*x/c)]; %#ok<SAGROW>
end

mode = C'*si;

for i = 1:5
    figure(i)
    plot(x, mode(i,:))
    hold on
end

main_mode = mode(1,:) + mode(2,:) + mode(3,:) + mode(4,:) + mode(5,:);

figure(6)
plot(x, main_mode);
hold on

%% orthogonality
%Ip=mean(Ipx);
Id=pi*(1.5)^4/32;
for i=1:5
    for j=1:5
        if i==j
            UO(i,j)=0;
        else
            UO(i,j)=(Ipx.*mode(i,:)*mode(j,:)')+Id*(mode(i,end)*mode(j,end));
        end
    end
end
disp("orthogonality");
disp(UO);
%% force
t=0:0.1:40;

a=sin(100*t);
b1=sin(wn(1)*t);
b2=sin(wn(2)*t);
b3=sin(wn(3)*t);
b4=sin(wn(4)*t);
b5=sin(wn(5)*t);

Zt1=median(diff(t))*conv(a,b1,'same');
Zt2=median(diff(t))*conv(a,b2,'same');
Zt3=median(diff(t))*conv(a,b3,'same');
Zt4=median(diff(t))*conv(a,b4,'same');
Zt5=median(diff(t))*conv(a,b5,'same');

Tur=1000;

Vx1=(Tur./Ipx).*(mode(1,:)./wn(1))*mode(1,401);
Vx2=(Tur./Ipx).*(mode(2,:)./wn(2))*mode(2,401);
Vx3=(Tur./Ipx).*(mode(3,:)./wn(3))*mode(3,401);
Vx4=(Tur./Ipx).*(mode(4,:)./wn(4))*mode(4,401);
Vx5=(Tur./Ipx).*(mode(5,:)./wn(5))*mode(5,401);


TETA1=Vx1'.*Zt1;
TETA2=Vx2'.*Zt2;
TETA3=Vx3'.*Zt3;
TETA4=Vx4'.*Zt4;
TETA5=Vx5'.*Zt5;
TETA=TETA1+TETA2+TETA3+TETA4+TETA5;


% Create a 3D surface plot
figure;
surf(t, x, TETA, 'EdgeColor', 'none');
xlabel('Time');
ylabel('X');
zlabel('Value');
title('3D Surface Plot: X vs. Time vs. TETA');
colormap('jet'); % You can choose a different colormap if needed
colorbar; % Add colorbar for reference





