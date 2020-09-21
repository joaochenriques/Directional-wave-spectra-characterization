clear
load('Condor_Clusters.mat')
Hm0=Results_Clusters.Hm0;
Te=Results_Clusters.Te;
Tp=Results_Clusters.Tp;
MeanDir=Results_Clusters.MeanDir;
energy=Results_Clusters.Power_spec;
time=Results_Clusters.time;
serie=datevec(time);
spring=find(serie(:,2)>=3 & serie(:,2) <=5);
summer=find(serie(:,2)>=6 & serie(:,2) <=8);
fall=find(serie(:,2)>=9 & serie(:,2) <=11);
winter=find(serie(:,2)==12 | serie(:,2) <=2);
density=1.025;
grav=9.81;

h=figure
a=1;
for i=1:12
    mes=find(serie(:,2)==i);
    mensal(mes)=i;
    maximo(i,:)=max(energy(mes));
    minimo(i,:)=min(energy(mes));
    media(i,:)=mean(energy(mes));
    mediana(i,:)=median(energy(mes));
    order=sort(energy(mes));
    cumulativo=cumsum(order);
    
    frequencia=linspace(0, max(energy(mes)),1000);
for j=0:length(frequencia)-1
    if j~0;
        frequencia(2,j)=sum(energy(mes)>frequencia(1,j) & energy(mes)<=frequencia(1,j+1));
    else
    end
end
    frequencia=(frequencia(2,:)/(sum(frequencia(2,:))))*100;
    frequencia=frequencia';
    fr_acum=linspace(0,max(energy(mes)),1000);
    fr_acum=fr_acum';
    fr_acum(:,2)=cumsum(frequencia);
    i25(1,:)=0:ceil(max(energy(mes)))/1000:ceil(max(energy(mes)));
    i25(2,:)=25;
    [q25(i,:),~]=curveintersect(i25(1,:),i25(2,:),fr_acum(:,1),fr_acum(:,2));
    i75(1,:)=0:ceil(max(energy(mes)))/1000:ceil(max(energy(mes)));
    i75(2,:)=75;
    [q75(i,:),~]=curveintersect(i75(1,:),i75(2,:),fr_acum(:,1),fr_acum(:,2));
    plot([a a],[q75(i,:) maximo(i,:)],'k');hold on
    plot([a a],[q25(i,:) minimo(i,:)],'k');
    plot([a-.3 a-.3 a+.3 a+.3 a-.3],[q25(i,:) q75(i,:) q75(i,:) q25(i,:) q25(i,:)],'k')
    plot([a-.4 a+.4],[mediana(i,:) mediana(i,:)],'k')
    plot(a,media(i,:),'ok')
    a=a+1;    
end
set(gca,'yscale','log','ytick',[0.1 1 5 10 20 50 100 200 500 1000])
set(gca,'ytick',[0 1 5 10 20 50 100 200 500 1000])
set(gca,'xtick',1:12, 'xticklabel',datestr(1:31:365,'mmm'))
grid on
ylim([0.1 2000])
xlim([0 13])
ylabel('P (kW/m)')
title('Energy Flux Seasonality')
set(gcf,'paperposition',[0 0 15 9]);
print(h,'-dpng','-r300','E_Sazonal.png')
close(h)

h=figure
a=1;
anos=2009:2017;
for i=1:9
    mes=find(serie(:,1)==anos(i));
    anual(mes)=i;
    maximo(i,:)=max(energy(mes));
    minimo(i,:)=min(energy(mes));
    media(i,:)=mean(energy(mes));
    mediana(i,:)=median(energy(mes));
    order=sort(energy(mes));
    cumulativo=cumsum(order);
    
    frequencia=linspace(0, max(energy(mes)),1000);
for j=0:length(frequencia)-1
    if j~0;
        frequencia(2,j)=sum(energy(mes)>frequencia(1,j) & energy(mes)<=frequencia(1,j+1));
    else
    end
end
    frequencia=(frequencia(2,:)/(sum(frequencia(2,:))))*100;
    frequencia=frequencia';
    fr_acum=linspace(0,max(energy(mes)),1000);
    fr_acum=fr_acum';
    fr_acum(:,2)=cumsum(frequencia);
    i25(1,:)=0:ceil(max(energy(mes)))/1000:ceil(max(energy(mes)));
    i25(2,:)=25;
    [q25(i,:),~]=curveintersect(i25(1,:),i25(2,:),fr_acum(:,1),fr_acum(:,2));
    i75(1,:)=0:ceil(max(energy(mes)))/1000:ceil(max(energy(mes)));
    i75(2,:)=75;
    [q75(i,:),~]=curveintersect(i75(1,:),i75(2,:),fr_acum(:,1),fr_acum(:,2));
    plot([a a],[q75(i,:) maximo(i,:)],'k');hold on
    plot([a a],[q25(i,:) minimo(i,:)],'k');
    plot([a-.3 a-.3 a+.3 a+.3 a-.3],[q25(i,:) q75(i,:) q75(i,:) q25(i,:) q25(i,:)],'k')
    plot([a-.4 a+.4],[mediana(i,:) mediana(i,:)],'k')
    plot(a,media(i,:),'ok')
    a=a+1;    
end
set(gca,'yscale','log','ytick',[0.1 1 5 10 20 50 100 200 500 1000])
set(gca,'xtick',1:9,'xticklabel',datestr(time(9000):365:time(end),'YYYY'))
grid on
ylim([0.1 2000])
xlim([0 10])
ylabel('P (kW/m)')
title('Energy Flux Interanual Variability')
set(gcf,'paperposition',[0 0 15 9]);
print(h,'-dpng','-r300','E_Interanual.png')
close(h)


distacumulada_generico(energy(mes),'lcb','Max. Energy Flux (kW/m)','titulo','Energy Flux Cumulative Distribution','variavel','P (kW/m)','printf',600,'sufix','_Energy_all')
distacumulada_generico(energy(winter),'lcb','Max. Energy Flux (kW/m)','titulo','Winter','variavel','P (kW/m)','printf',600,'sufix','_Energy_winter')
distacumulada_generico(energy(summer),'lcb','Max. Energy Flux (kW/m)','titulo','Summer','variavel','P (kW/m)','printf',600,'sufix','_Energy_summer')
distacumulada_generico(energy(fall),'lcb','Max. Energy Flux (kW/m)','titulo','Fall','variavel','P (kW/m)','printf',600,'sufix','_Energy_fall')
distacumulada_generico(energy(spring),'lcb','Max. Energy Flux (kW/m)','titulo','Spring','variavel','P (kW/m)','printf',600,'sufix','_Energy_spring')


%Distribuições conjuntas
ccc=flipud(hot);
ccc=ccc(14:end,:);

diag_direc_OTP(Te,MeanDir,'dirl',4,'bt',2,'nclass',50,'ps',3,'un','s','eixo',[4 8 12 16],'el',0.5,'titulo','MeanDir X Te','printf',300,'sufix','_TeDir_paper_hot','lcb','Occurrence (%)','ccc',ccc);
diag_direc_OTP(Hm0,MeanDir,'dirl',4,'bt',2,'nclass',50,'ps',3,'un','m','eixo',[2 5 8 11],'el',0.5,'titulo','MeanDir X Hm0','printf',300,'sufix','_Hm0Dir_paper_hot','lcb','Occurrence (%)','ccc',ccc);


classe_h=0.1:.2:max(Hm0)+4;
classe_t=2:0.2:max(Te)+4;
for i=1:length(classe_h)-1
    h_classe=find(Hm0>classe_h(i) & Hm0<classe_h(i+1));
    Hm(i)=mean(Hm0(h_classe));
    for j=1:length(classe_t)-1
        Hi(i)=(classe_h(i) + classe_h(i+1))/2;
        Ti(j)=(classe_t(j) + classe_t(j+1))/2;
        ind=find(Te(h_classe)>classe_t(j) & Te(h_classe)<classe_t(j+1));
        prob(i,j)=length(find(Te(h_classe)>classe_t(j) & Te(h_classe)<classe_t(j+1)))/length(Te);
        E_class(i,j)=(1/(64*pi))*density*grav^2*Hi(i)^2*Ti(j);
    end
end
prob(prob==0)=nan;
prob=100.*prob;

h=figure
pcolor(classe_t(1:end-1),classe_h(1:end-1),prob);shading flat; grid on; hold on
contour(classe_t(1:end-1),classe_h(1:end-1),E_class,[5 20 50 100 200 400 800],'k')
yy=colorbar;colormap(ccc);
title(yy,'Occur.(%)','fontsize',10)
ylabel('Hm0 (m)','fontsize',10)
xlabel('Te (s)','fontsize',10)
title('Te X Hm0','fontsize',12)
xlim([2 16])
ylim([0 12])
caxis([0 0.5])
text(2.2,1.3,'5 kW/m');
text(2.9,2.6,'20 kW/m');
text(4,3.8,'50 kW/m');
text(5,5.1,'100 kW/m');
text(6.4,6.7,'200 kW/m');
text(8.1,8.7,'400 kW/m');
text(10.4,11.25,'800 kW/m');
set(gcf,'paperposition',[0 0 18 10]);
print(h,'-dpng','-r300','Hm0_Tp_paper_hot.png')
close(h)


%% Figura Exemplo Spectro-Parametro

load('Spectra_Condor.mat');
T=1./f;
time=[time time(end)+1/48:1/48:time(end)+1-1/48];
load('Condor_Clusters.mat')
Hm0=Results_Clusters.Hm0;
Tp=Results_Clusters.Tp;
MeanDir=Results_Clusters.MeanDir;
energy=Results_Clusters.Power_spec;
ccc=flipud(hot);
ccc=ccc(14:end,:);


sample=7578;
espectro=spec(:,:,sample);
% Ajustando espectro idealizado ao exemplo
SB = bretschneider(f*2*pi,[Hm0(sample) Tp(sample)]);
D = spreading(24,'cos2s',deg2rad(MeanDir(sample)),[],SB.w,0);
SBd = mkdspec(SB,D); % Directional spectrum
SBd.S(SBd.S==0)=nan;

theta_exemplo=SBd.theta;
theta_exemplo(theta_exemplo<0)=theta_exemplo(theta_exemplo<0)+(2*pi);
theta_exemplo=[theta_exemplo(13:end);theta_exemplo(1:12)];
S_exemplo=[SBd.S(13:end,:);SBd.S(1:12,:)]*2*pi;

h=figure
polar_spect_part2(T,theta,espectro);colorbar;colormap(flipud(hot));shading flat;hold on
[c,g]=polarcont2((2*pi)./SBd.w,deg2rad(cart2naut(rad2deg(theta_exemplo))),S_exemplo',8,'r');
% clabel(c,h)
hcb=colorbar;
s = sprintf('m^2/s/%c', char(176));
title(hcb,{'E(f,\theta)',s})
lim=[4 8 12 16 20];
theta2=[theta;theta(1)];
for g=1:5
    for j=1:length(theta)
plot([lim(g)*sin(theta2(j)) lim(g)*sin(theta2(j+1))],[lim(g)*cos(theta2(j)) lim(g)*cos(theta2(j+1))],':k');hold on
    end
end
text(.7,3,'4s','fontsize',8,'color','k');
text(.7,7.1,'8s','fontsize',8,'color','k');
text(.7,11.1,'12s','fontsize',8,'color','k');
text(.7,15.1,'16s','fontsize',8,'color','k');
text(.7,19.1,'20s','fontsize',8,'color','k');
xlim([-21 21])
ylim([-21 21])
title('Energy Density Spectra - 05/12/2018 20:30')
axis off
set(gcf,'paperposition',[0 0 13 13]);
print(h,'-dpng','-r300','Exemplo_espectro2.png')
close(h)

h=figure
subplot(1,2,1)
plot(T,trapz(theta_exemplo,espectro'),'k','linewidth',2);hold on
plot(T,trapz(theta_exemplo,S_exemplo),':r','linewidth',2)
grid on;xlim([2 25])
xlabel('Period (s)')
ylabel('E(f) (m^2/s)')
legend('Real Spectra','Idealized Spectra','location','northoutside')
subplot(1,2,2)
plot(rad2deg(theta_exemplo),trapz(f,espectro),'k','linewidth',2);hold on
plot(rad2deg(theta_exemplo),trapz(f,S_exemplo'),':r','linewidth',2)
grid on;xlim([0 360])
xlabel('Direction (º)')
ylabel('E(\theta) (m^2/º)')
legend('Real Spectra','Idealized Spectra','location','northoutside')
set(gcf,'paperposition',[0 0 16 8]);
print(h,'-dpng','-r300','Exemplo_espectro_1D.png')
close(h)


alfa=(5*Hm0(sample)^2*Tp(sample)^-4)/(16*(9.81^2)*(2*pi)^-4)

for i=1:length(f)
        PM_f(i)=alfa*(9.81^2)*((2*pi)^-4)*(f(i)^-5)*exp(-(5/4)*((f(i)*Tp(sample))^-4));
end

ww=T/(2*pi);
S_interp=interp1(ww,PM_f,SB.w);
figure
plot(T,SB.S*2*pi,'r',T,PM_f,':k')


% Calcula porcentagem de cada cluster nas estacoes do ano
for i=1:9
    OCC(i,:)=round((length(find(Results_Clusters.idx==i))/length(Results_Clusters.idx))*100,2);
    Oc_winter(i,:)=round((length(find(Results_Clusters.idx(winter)==i))/length(winter))*100,2);
    Oc_spring(i,:)=round((length(find(Results_Clusters.idx(spring)==i))/length(spring))*100,2);
    Oc_summer(i,:)=round((length(find(Results_Clusters.idx(summer)==i))/length(summer))*100,2);
    Oc_fall(i,:)=round((length(find(Results_Clusters.idx(fall)==i))/length(fall))*100,2);
end
rr=[5 1 8 4 3 7 6 2 9];
[~,order]=sort(rr);
tabela=[rr; OCC';Oc_winter';Oc_spring';Oc_summer';Oc_fall'];
tabela=tabela(:,order);
    


% Calcula porcentagem de cada cluster meses e anos

timeserie=datevec(Results_Clusters.time);
anos=2009:2017;
meses=1:12;
cc=[1 0 0;.7 0 0;.4 0 0;0 0 1;0 0 .7;0 0 .4;0 1 0;.0 .7 0;0 .4 0];
ff=cc([1 4 5]);
for i=1:12 
    [a,b]=find(timeserie(:,2)==meses(i));
    idx_mes=Results_Clusters.idx(a,:);
    power_spec_mes=Results_Clusters.Power_param(a,:);
    for j=1:9
        [e,~]=find(idx_mes==j);
        ocorr_mes(j,i)=length(e)./length(idx_mes);
    end
    [u,months(i,:)]=month(datenum([2010 i 1 1 1 1]));
end

for i=1:9 
    [a,b]=find(timeserie(:,1)==anos(i));
    idx_ano=Results_Clusters.idx(a,:);
    power_spec_ano=Results_Clusters.Power_param(a,:);
    for j=1:9
        [e,~]=find(idx_ano==j);
        ocorr_ano(j,i)=length(e)./length(idx_ano);
    end
    E_real_anos(i,:)=mean(power_spec_ano);
end


 % reordenando espectros para ficarem associados pelo tipo no gráfico
 order2=[4 9 6 1 2 8 5 3 9];
ocorr_mes=ocorr_mes(order2,:);
ocorr_ano=ocorr_ano(order2,:);

h=figure
area(ocorr_mes','linestyle','none');colormap(cc);
ylim([0 1]);xlim([1 12])
hcb=legend('1','2','3','4','5','6','7','8','9','Orientation','horizontal','location','southoutside');
title('Representative Spectra #');
legend('boxoff')
title('Spectra # Occurrence by Month (%)')
set(gca,'ytick',[])
set(gca,'xticklabel',months);

set(gcf,'paperposition',[0 0 20 14]);
print(h,'-dpng','-r300','Porcentagem_clusters_S.png')
close(h)

cent=Results_Clusters.Centroids_Originais;
for i=1:9
    [e,r]=find(Results_Clusters.idx==i);
    hs(i)=4*sqrt(trapz(Results_Clusters.theta,trapz(Results_Clusters.f,mean(Results_Clusters.Spectra(:,:,e),3))));
    ocorr(i,:)=length(e)./length(Results_Clusters.idx); 
    h=figure
    polarcont_naut(Results_Clusters.T,Results_Clusters.theta,mean(Results_Clusters.Spectra(:,:,e),3));colorbar;colormap(flipud(hot));shading flat
end
% Plot dos centroids originais e espectros idealizados para cada um deles

cent=Results_Clusters.Centroids_Originais;

for i=1:9
    hs(i)=4*sqrt(trapz(theta,trapz(f,cent(:,:,i))));
    f_spec(i,:)=trapz(theta,cent(:,:,i)');
[~,b]=max(f_spec(i,:));
if b<23
p = polyfit(f(b-2:b+2)',f_spec(i,b-2:b+2),2);
x=f(b-2):0.001:f(b+2);
elseif b<24
p = polyfit(f(b-2:b+1)',f_spec(i,b-2:b+1),2);
x=f(b-2):0.001:f(b+1);
else
p = polyfit(f(b-2:b)',f_spec(i,b-2:b),2);
x=f(b-2):0.001:f(b);
end
[a,b]=max(polyval(p,x));
Tp(i,:)=1/x(b);
SB = bretschneider(f*2*pi,[hs(i) Tp(i,:)]);
h=figure
plot(T,f_spec(i,:),'k');hold on
plot(T,SB.S.*2*pi,'r');grid on
title(num2str(hs(i)))
end


    
    