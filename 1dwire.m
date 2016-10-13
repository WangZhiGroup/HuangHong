clear all
tic
rand('state',sum(100*clock));
%###############################
tms=300;
%Parameter Initialization
nsize=1000;layers=1;Sites=nsize/2;
n=4*nsize*layers;
dphi=0.2*pi/tms;
gapa=0.5;couple=0.005;t=11;a=0.67; mu=0.2;h_t=couple*t;
v=linspace(0,1,tms+1);
espectrum1=zeros(n,tms+1);
espectrum2=zeros(n,tms+1);

for k=1:tms+1
    disp(k)
    V=[v(k),0,0];
    Hc1a=Hamiltonian_NW(nsize/2,layers,t,mu,a,V,gapa*exp(1i*pi/2));
    Hc1b=Hamiltonian_NW(nsize/2,layers,t,mu,a,V,gapa*exp(1i*(pi/2+dphi)));
    Hc2=Hamiltonian_NW(nsize/2,layers,t,mu,a,V,gapa);
    
    Hca = [Hc1a,zeros(4*Sites*layers);zeros(4*Sites*layers),Hc2];
    Hca(4*layers*Sites - 4*layers+1:4*Sites*layers,4*layers*Sites + 1:4*layers*Sites + 4*layers) = h_t*kron(eye(layers),diag([-1,-1,1,1]));
    Hca(4*layers*Sites + 1:4*layers*Sites + 4*layers,4*layers*Sites - 4*layers+1:4*layers*Sites) = h_t*kron(eye(layers),diag([-1,-1,1,1]));
    
    Hcb = [Hc1b,zeros(4*Sites*layers);zeros(4*Sites*layers),Hc2];
    Hcb(4*layers*Sites - 4*layers+1:4*Sites*layers,4*layers*Sites + 1:4*layers*Sites + 4*layers) = h_t*kron(eye(layers),diag([-1,-1,1,1]));
    Hcb(4*layers*Sites + 1:4*layers*Sites + 4*layers,4*layers*Sites - 4*layers+1:4*layers*Sites) = h_t*kron(eye(layers),diag([-1,-1,1,1]));
    
    [~,evalue1]=eig(Hca);
    [~,evalue2]=eig(Hcb);
    
    espectrum1(:,k)=diag(evalue1);
    espectrum2(:,k)=diag(evalue2);
    
    
end

jall=sum((espectrum1(n/2+1:n,:)-espectrum2(n/2+1:n,:))./dphi);

plot(v,espectrum1,'b')
title('Energy Spectrum1(V=1.0)','FontSize',14)
xlabel('\theta/\pi');
ylabel('Energy');
set(gca,'xlim',[0,2])
set(gca,'ylim',[-0.5,0.5])
saveas(gcf,'esp1.eps','psc2')

%  plot(hz,espectrum2,'b')
%  title('Energy Spectrum2(V=1.0)','FontSize',14)
%  xlabel('V_x');
%  ylabel('Energy');
%  set(gca,'xlim',[0,2])
%  set(gca,'ylim',[-0.5,0.5])
%  saveas(gcf,'esp2.eps','psc2')

plot(v,jall,'LineWidth',2)
xlabel('V_x');
ylabel('Current');
set(gca,'xlim',[0,2])
grid
saveas(gcf,'jall.eps','psc2')

save('espectrum1.mat','espectrum1')
save('espectrum2.mat','espectrum2')
save('jall.mat','jall')


clear evector1 evalue1;
clear hmtd hmtv1 hmto a1;

toc

