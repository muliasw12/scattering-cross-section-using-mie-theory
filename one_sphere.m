clearvars;
clc;
tic
%program ini memakai rumusan teori Mie untuk satu bola, konvergen di N = 2

%parameter-parameter 
E0 = 1;
c = 3E+08;
mu0 = 4*pi*1E-07;
eps0 = 8.85*1E-12;
R = 100; %in nanometer
N = 1; %medium
Nt = 4; %number of Bessel mode
% load J_C
% lambda = 1243./J_C(:,1); %in nanometer
load Si
lambda = 1000*Si(:,1); %in nanometer
lambdamin = 400; lambdamax = 1000;
dex = find(lambda>lambdamin & lambda<lambdamax);
Lambda = lambda(min(dex):max(dex));
lambda = lambdamin:0.5:lambdamax; %in nanometer
f = c./(lambda*1E+3); %frequency in THz
k = 2*pi*N./(lambda);
x = k*R; %size parameter
% m = (J_C(:,6)+(1i*J_C(:,7)))/N; %for gold
%m = (J_C(:,4)+(1i*J_C(:,5)))/N; %for silver
%m = (J_C(:,2)+(1i*J_C(:,3)))/N; %for copper
m = (Si(:,2)+(1i*Si(:,3)))/N; %for silicon
mu1 = 1;
% m = (2.41+0*1i*0.030209)*ones(length(lambda),1)/N; %for dielectrics
M = m(min(dex):max(dex));
m = interp1(Lambda,M,lambda,'pchip');

% m = (4+1i*0.03)*ones(1,size(m,2)); % Pak Alex

Lambda = 694; %ingin meninjau pada lambda berapa
in = find(lambda==Lambda); 

%permittivity from Drude model
% epsi = 9; hwp = 9; hgamma = 0.05; j=1;
% for i=lambdamin:1:lambdamax
%     hw = 1243/i;
%     eps = epsi - hwp^2/(hw*(hw+1i*hgamma));
%     m(j) = sqrt(eps)/N;
%     j = j+1;
% end

%fungsi-fungsi yg dipakai
j = @(n,x) sqrt(pi./(2*x)).*besselj(n+0.5,x); %bessel sferis
h = @(n,x) sqrt(pi./(2*x))*(besselj(n+0.5,x) + 1i*bessely(n+0.5,x)); %hankel sferis
psi = @(n,x) x.*j(n,x); %riccati-bessel dan turunannya
dpsi = @(n,x) j(n,x) + (x/(2*n+1))*(n*j(n-1,x) - (n+1)*j(n+1,x));
xi = @(n,x) x.*h(n,x);
dxi = @(n,x) h(n,x) + (x/(2*n+1))*(n*h(n-1,x) - (n+1)*h(n+1,x));
%hitung koefisien an dan bn
for i=1:length(lambda)
    for n=1:Nt
        a(n,i) = (m(i)*psi(n,m(i)*x(i))*dpsi(n,x(i))-mu1*psi(n,x(i))*dpsi(n,m(i)*x(i)))...
            /(m(i)*psi(n,m(i)*x(i))*dxi(n,x(i))-mu1*xi(n,x(i))*dpsi(n,m(i)*x(i)));
        b(n,i) = (mu1*psi(n,m(i)*x(i))*dpsi(n,x(i))-m(i)*psi(n,x(i))*dpsi(n,m(i)*x(i)))...
            /(mu1*psi(n,m(i)*x(i))*dxi(n,x(i))-m(i)*xi(n,x(i))*dpsi(n,m(i)*x(i)));
    end
    s1 = 0; s2 = 0;
    for n=1:Nt
        s1 = s1 + (2*n+1)*((abs(a(n,i)))^2 + (abs(b(n,i)))^2);
        s2 = s2 + (2*n+1)*real(a(n,i) + b(n,i));
    end
    Csca(i) = (2*pi/(k(i)^2))*s1; Qsca(i) = Csca(i)/(pi*R^2);
    Cext(i) = (2*pi/(k(i)^2))*s2; Qext(i) = Cext(i)/(pi*R^2);
    Cabs(i) = Cext(i) - Csca(i); Qabs(i) = Cabs(i)/(pi*R^2);
    sde = (2*1+1)*(abs(a(1,i)))^2; %dipol e
    sdm = (2*1+1)*(abs(b(1,i)))^2; %dipol m
    Csde(i) = (2*pi/(k(i)^2))*sde; Qsde(i) = Csde(i)/(pi*R^2);
    Csdm(i) = (2*pi/(k(i)^2))*sdm; Qsdm(i) = Csdm(i)/(pi*R^2);  
end
figure()
plot(lambda,Csca*1E-6,'Color','k','LineWidth',1)
xlabel('wavelength (nm)'); ylabel('C_{sca} (um^2)')
title (['silicon, r = ', num2str(R), ' nm, dispersif'])

figure()
plot(lambda,Qsca,'Color','k','LineWidth',1)
xlabel('wavelength (nm)'); ylabel('Q_{sca}')
title (['silicon, r = ', num2str(R), ' nm, dispersif'])

% figure
% plot(lambda,Csca,'Color','k','LineWidth',2)
% xlabel('\lambda (nm)'); ylabel('cross section (nm^2)')
% title (['gold nanosphere ', num2str(R), ' nm'])
% hold on
% plot(lambda,Cext,'Color','g','LineWidth',2)
% hold on
% plot(lambda,Cabs,'Color','r','LineWidth',2)
% legend ('C_{sca}','C_{ext}','C_{abs}')

% figure
% plot(lambda,Qsca,'Color','k','LineWidth',2)
% xlabel('\lambda (nm)'); ylabel('efficiency')
% title (['gold nanosphere ', num2str(R), ' nm'])
% hold on
% plot(lambda,Qext,'Color','g','LineWidth',2)
% hold on
% plot(lambda,Qabs,'Color','r','LineWidth',2)
% legend ('Q_{sca}','Q_{ext}','Q_{abs}')

% plot(f,Csca,'Color','k','LineWidth',2)
% xlabel('f (THz)'); ylabel('cross section (nm^2)')
% title (['gold nanosphere ', num2str(R), ' nm'])
% hold on
% plot(f,Cext,'Color','g','LineWidth',2)
% hold on
% plot(f,Cabs,'Color','r','LineWidth',2)
% legend ('C_{sca}','C_{ext}','C_{abs}')

% %plot intensitas medan listrik
% %bangun grid
% w = 100; del = 1; %ukuran grid
% [x,y] = meshgrid(-w:del:w); z = 0*ones(size(x));
% r = sqrt(x.^2+y.^2+z.^2);
% the = acos(z./r);
% phi = atan2(y,x);
% for i=1:2*w/del+1
%     for j=1:2*w/del+1
%         if r(i,j) >= R %5/k(in) %plot medan terhambur di luar bola
%             for n=1:Nt
%                 %vektor M dan N, r = radial, t = theta, p = phi
%                 Me1nt = -sin(phi(i,j))*pin(n,cos(the(i,j)))*h(n,k(in)*r(i,j));
%                 Me1np = -cos(phi(i,j))*taun(n,cos(the(i,j)))*h(n,k(in)*r(i,j));
%                 Mo1nt = cos(phi(i,j))*pin(n,cos(the(i,j)))*h(n,k(in)*r(i,j));
%                 Mo1np = -sin(phi(i,j))*taun(n,cos(the(i,j)))*h(n,k(in)*r(i,j));
%                 Ne1nr = h(n,k(in)*r(i,j))*cos(phi(i,j))*n*(n+1)*sin(the(i,j))*pin(n,cos(the(i,j)))/(k(in)*r(i,j));
%                 Ne1nt = cos(phi(i,j))*taun(n,cos(the(i,j)))*dxi(n,k(in)*r(i,j))/(k(in)*r(i,j));
%                 Ne1np = -sin(phi(i,j))*pin(n,cos(the(i,j)))*dxi(n,k(in)*r(i,j))/(k(in)*r(i,j));
%                 No1nr = h(n,k(in)*r(i,j))*sin(phi(i,j))*n*(n+1)*sin(the(i,j))*pin(n,cos(the(i,j)))/(k(in)*r(i,j));
%                 No1nt = sin(phi(i,j))*taun(n,cos(the(i,j)))*dxi(n,k(in)*r(i,j))/(k(in)*r(i,j));
%                 No1np = cos(phi(i,j))*pin(n,cos(the(i,j)))*dxi(n,k(in)*r(i,j))/(k(in)*r(i,j));
%                 %medan terhambur E dan H
%                 En = (1i)^n*E0*(2*n+1)/(n*(n+1));
%                 eta = N/(c*mu0);
%                 Ern(n) = En*(1i)*a(n,in)*Ne1nr; %Hrn(n) = En*(1i)*b(n,in)*No1nr;
%                 Etn(n) = En*(-b(n,in)*Mo1nt + (1i)*a(n,in)*Ne1nt); %Htn(n) = En*(a(n,in)*Me1nt + (1i)*b(n,in)*No1nt);
%                 Epn(n) = En*(-b(n,in)*Mo1np + (1i)*a(n,in)*Ne1np); %Hpn(n) = En*(a(n,in)*Me1np + (1i)*b(n,in)*No1np);         
%             end
%             Er(i,j) = sum(Ern); Et(i,j) = sum(Etn); Ep(i,j) = sum(Epn);
%             %Hr(i,j) = eta*sum(Hrn); Ht(i,j) = eta*sum(Htn); Hp(i,j) = eta*sum(Hpn);
%             Ex(i,j) = Er(i,j)*sin(the(i,j))*cos(phi(i,j)) + Et(i,j)*cos(the(i,j))*cos(phi(i,j)) - Ep(i,j)*sin(phi(i,j));
%             Ey(i,j) = Er(i,j)*sin(the(i,j))*sin(phi(i,j)) + Et(i,j)*cos(the(i,j))*sin(phi(i,j)) + Ep(i,j)*cos(phi(i,j));
%             Ez(i,j) = Er(i,j)*cos(the(i,j)) - Et(i,j)*sin(the(i,j));
%             IE(i,j) = (abs(Er(i,j)))^2 + (abs(Et(i,j)))^2 + (abs(Ep(i,j)))^2; %intensitas
%             %IH(i,j) = (abs(Hr(i,j)))^2 + (abs(Ht(i,j)))^2 + (abs(Hp(i,j)))^2;
%         end
%     end
% end
% figure;
% surf(x,y,IE);colormap jet;shading interp;
% hold on
% view(0,90);colorbar;
% xlabel ('x (nm)');
% ylabel ('y (nm)');
% title ({['intensitas radiasi medan E terhambur pada bola'];['R = ', num2str(R), ' nm, \lambda = ', ...
%     num2str(Lambda), ' nm, I_0 = ', num2str(E0^2), ' W/m^2']})
% figure;
% quiver(x,z,Ex,Ez)
% 
% figure;
% surf(z,x,IH);colormap jet;shading interp;
% hold on
% view(0,90);colorbar;
% xlabel ('z (nm)');
% ylabel ('x (nm)');
% title ({['intensitas radiasi medan H terhambur pada bola'];['R = ', num2str(R), ' nm, \lambda = ', ...
%     num2str(Lambda), ' nm, I_0 = ', num2str(E0^2), ' W/m^2']})

% %polar plot
% theta = 0:pi/100:2*pi; phi_ = 0; r_ = 2000; th = 1;
% S1 = zeros(1,length(theta)); S2 = zeros(1,length(theta)); 
% for i=1:length(theta) %scattering amplitude
%     if theta(i) > pi & th == 1
%         th = th + 1;
%         phi_ = phi_ + pi;
%     end
%     for n=1:N
%         S1(i) = S1(i) + (2*n+1)*(a(n,in)*pin(n,cos(theta(i))) + b(n,in)*taun(n,cos(theta(i))))/(n*(n+1));
%         S2(i) = S2(i) + (2*n+1)*(a(n,in)*taun(n,cos(theta(i))) + b(n,in)*pin(n,cos(theta(i))))/(n*(n+1));
%     end
%     Est(i) = E0*exp(1i*k(in)*r_)*cos(phi_)*S2(i)/(-1i*k(in)*r_);
%     Esp(i) = -E0*exp(1i*k(in)*r_)*sin(phi_)*S1(i)/(-1i*k(in)*r_);
%     %I1(i) = (abs(S1(i)))^2; %intensity
%     %I2(i) = (abs(S2(i)))^2;
%     I(i) = (abs(Est(i)))^2 + (abs(Esp(i)))^2;
% end
% figure
% polarplot(theta,I)
% title (['intensitas radiasi pada \phi = ',num2str(phi_)])
toc

function y = pin(n,x)
if x == 1
    x = x - 1E-05;
elseif x == -1
    x = x + 1E-05;
end
y = P(1,n,x)/sqrt(1-x^2);
end

function y = taun(n,x)
if x == 1
    x = x - 1E-05;
elseif x == -1
    x = x + 1E-05;
end
y = dP(1,n,x);
end