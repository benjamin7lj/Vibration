clc
clear
close all;
%% data
G=8*10^10; % shear modulus
p=7800;% density 
D1=0.5; % the diameter of beam is d0+d1(1-x/l)
D0=0.25;% the diameter of beam is d0+d1(1-x/l)
L=4;% l= 4 (length of beam)
Kt=5*10^6;%k=50Mnm
x=0:0.01:L;
D=D0+D1*(1-x/L); % diameter D of beam
Jx=pi*D.^4/32;
Dg=1.5; % diameter of gear is 60 cm
Jg=p*pi*Dg^4/32; 
%% galerkin(w)
c=sqrt(G/p);
syms w
f=G*Jx(end)*(w)*cos(w)/(Jg*w^2-Kt);
g=sin(w);
modes = zeros(1, 5);
for mode = 1:5
    initial=1 +pi* (mode - 1);

    i = vpasolve(f == g, w, initial);
    if i>0
      modes(mode) = i;
    end
end
modes=sort(modes);
disp('phi_s "galerkin" for natural boundary condition');
disp(modes);