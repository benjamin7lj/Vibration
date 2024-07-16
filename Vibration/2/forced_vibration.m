%this is the main code for Forced vibration of membrane

clc
clear
close all;

syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12...
    c13 c14 c15 c16 x y w t THA
C=[c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 ...
    c11 c12 c13 c14 c15 c16];
a=0.1;b=0.1;
ro=2700;
p1=1;p2=-0.5;
p3=1;p4=-0.5;
p5=0.7;
Px=2000*(p1+p2*(y/b));
Py=2000*(p3+p4*(x/a));
Txy=p5*1000;
F0=10^8;
F01=1;
w_f=70;

W(1)=c1*sin((x*pi)/a)*sin((y*pi)/b);
W(2)=c2*sin((1*x*pi)/a)*sin((2*y*pi)/b);
W(3)=c3*sin((2*x*pi)/a)*sin((1*y*pi)/b);
W(4)=c4*sin((2*x*pi)/a)*sin((2*y*pi)/b);
W(5)=c5*sin((1*x*pi)/a)*sin((3*y*pi)/b);
W(6)=c6*sin((3*x*pi)/a)*sin((1*y*pi)/b);
W(7)=c7*sin((2*x*pi)/a)*sin((3*y*pi)/b);
W(8)=c8*sin((3*x*pi)/a)*sin((2*y*pi)/b);
W(9)=c9*sin((1*x*pi)/a)*sin((4*y*pi)/b);
W(10)=c10*sin((4*x*pi)/a)*sin((y*pi)/b);
W(11)=c11*sin((2*x*pi)/a)*sin((4*y*pi)/b);
W(12)=c12*sin((4*x*pi)/a)*sin((2*y*pi)/b);
W(13)=c13*sin((3*x*pi)/a)*sin((3*y*pi)/b);
W(14)=c14*sin((4*x*pi)/a)*sin((3*y*pi)/b);
W(15)=c15*sin((3*x*pi)/a)*sin((4*y*pi)/b);
W(16)=c16*sin((4*x*pi)/a)*sin((4*y*pi)/b);

W1(1)=sin((x*pi)/a)*sin((y*pi)/b);
W1(2)=sin((1*x*pi)/a)*sin((2*y*pi)/b);
W1(3)=sin((2*x*pi)/a)*sin((1*y*pi)/b);
W1(4)=sin((2*x*pi)/a)*sin((2*y*pi)/b);
W1(5)=sin((1*x*pi)/a)*sin((3*y*pi)/b);
W1(6)=sin((3*x*pi)/a)*sin((1*y*pi)/b);
W1(7)=sin((2*x*pi)/a)*sin((3*y*pi)/b);
W1(8)=sin((3*x*pi)/a)*sin((2*y*pi)/b);
W1(9)=sin((1*x*pi)/a)*sin((4*y*pi)/b);
W1(10)=sin((4*x*pi)/a)*sin((y*pi)/b);
W1(11)=sin((2*x*pi)/a)*sin((4*y*pi)/b);
W1(12)=sin((4*x*pi)/a)*sin((2*y*pi)/b);
W1(13)=sin((3*x*pi)/a)*sin((3*y*pi)/b);
W1(14)=sin((4*x*pi)/a)*sin((3*y*pi)/b);
W1(15)=sin((3*x*pi)/a)*sin((4*y*pi)/b);
W1(16)=sin((4*x*pi)/a)*sin((4*y*pi)/b);


phi=sum(W);

dif_x_phi=diff(phi,x);
dif_y_phi=diff(phi,y);

F_U=Px*(dif_x_phi)^2+Py*(dif_y_phi)^2+...
    2*Txy*(dif_x_phi)*(dif_y_phi);
F_U_int_x=int(F_U,x,0,a);
U_max=0.5*(int(F_U_int_x,y,0,b));

F_T=ro*phi^2;
F_T_int_x=int(F_T,x,0,a);
T_max=0.5*int(F_T_int_x,y,0,b);

M = vpa(zeros(16,16));
K = vpa(zeros(16,16));
for i=1:16
    U_temp = diff(U_max,C(i));
    T_temp = diff(T_max,C(i));
    for j = 1:16
        K(i,j) = diff(U_temp,C(j));
        M(i,j) = diff(T_temp,C(j));
    end
end


[v,vw]=eig(double(K),double(M));
w_Ritz=(diag(vw)).^0.5;
for i=1:16
    v(:,i) = v(:,i)/(max(abs(v(:,i))));
end

for i=1:16
    PHI1(i)=W1*v(:,i);
    figure
%     fsurf(PHI1(i),[0,1])
%     hold on
    fcontour(PHI1(i),[0,a],'Fill','on')
    hold on
end

%% orthogonality
ORT=zeros(16,16);

for i=1:16
    for j=1:16
        ORT(i,j)=v(:,i)'*M*v(:,j);
    end
end  

%% forced
% F=-0.5*F01*sin(w_f*t)*dirac(x-(a/2))*dirac(y-(b/2))+...
%     F01*sin(w_f*t)*dirac(x-(a/4))*dirac(y-(b/4))+...
%     F01*sin(w_f*t)*dirac(x-(3*a/4))*dirac(y-(3*b/4))+...
%     F01*sin(2*w_f*t)*dirac(x-(3*a/4))*dirac(y-(b/4))+...
%     F01*sin(2*w_f*t)*dirac(x-(a/4))*dirac(y-(3*b/4));
F=F0*sin(w_f*t)*(x^2-a*x)*(y^2-b*y);
for i =1:12
    NN=F*PHI1(i);
    int_x_NN=int(NN,x,0,a);
    N(i)=int(int_x_NN,y,0,b);
end

for i =1:12
    Nb(i)=subs(N(i),t,THA)*sin(w_Ritz(i)*(THA-t));
    NB(i)=int(Nb(i),THA,0,t);
    Wxyt(i)=(1/w_Ritz(i))*PHI1(i)*NB(i);
end

Wxyt_total=sum(Wxyt);
TTT=0:0.05:1.05;
for i=1:length(TTT)
    WWxyt(i)=subs(Wxyt_total,t,TTT(i));
    figure
    fsurf(WWxyt(i),[0,a])
    %hold on
    figure
    fcontour(WWxyt(i),[0,a],'Fill','on')
    pause(0.1)
end
