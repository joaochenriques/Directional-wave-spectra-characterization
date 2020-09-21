% Rotina para processar dados espectrais obtidos em 10 anos de WW3
% Por Pedro Ribeiro, em 30/08/2018

clear
close all
density=1.025;
grav=9.81;

S=load('Spectra_Condor.mat');

% Interpolando classes de direção
theta2=[S.theta;S.theta(1)+S.theta(end)];
med_T=S.theta(2:2:end);
for j=1:length(S.spec)
    spec_aux=[S.spec(:,:,j) S.spec(:,1,j)];
for i=1:length(S.theta)/2
    specT(:,i,j)=(trapz(theta2(2*i-1:2*i+1),spec_aux(:,2*i-1:2*i+1)'))./(theta2(2*i+1)-theta2(2*i-1));
end
    clear spec_aux
end

indices=[8 10 14 18];
            
% Dividindo o espectro em quatro classes de
% frequencia e Gerando "Hs equivalente"

for i=1:4
    T_aux=find(S.T<indices(i));
    T_ind(i)=T_aux(1);
end

Ti=[S.T(1:T_ind(4)-1,:);indices(4); S.T(T_ind(4):T_ind(3)-1,:);indices(3); S.T(T_ind(3):T_ind(2)-1,:);indices(2); S.T(T_ind(2):T_ind(1)-1,:);indices(1);S.T(T_ind(1):end,:)];
fi=[1./Ti(1:end-1);S.f(end)];
speci=ones(length(fi),length(med_T),length(specT))*nan;       
for i=1:length(specT)
    for j=1:length(med_T)
speci(:,j,i)=interp1(S.f,specT(:,j,i),fi);
    end
end


ind1=find(fi==1/indices(4));
ind2=find(fi==1/indices(3));
ind3=find(fi==1/indices(2));
ind4=find(fi==1/indices(1));
fff=[fi(1) fi(ind1) fi(ind2) fi(ind3) fi(ind4) fi(end)];
med=[mean([fi(1);fi(ind1)]) mean([fi(ind1);fi(ind2)]) mean([fi(ind2);fi(ind3)]) mean([fi(ind3);fi(ind4)]) mean([fi(ind4);fi(end)])];

for i=1:length(speci)
He_1=trapz(fi(1:ind1),speci(1:ind1,:,i)).*[.26 diff(med_T)'];
He_2=trapz(fi(ind1:ind2),speci(ind1:ind2,:,i)).*[.26 diff(med_T)'];
He_3=trapz(fi(ind2:ind3),speci(ind2:ind3,:,i)).*[.26 diff(med_T)'];
He_4=trapz(fi(ind3:ind4),speci(ind3:ind4,:,i)).*[.26 diff(med_T)'];
He_5=trapz(fi(ind4:end),speci(ind4:end,:,i)).*[.26 diff(med_T)'];
Hs(i,:)=4.*sqrt(trapz(med_T,trapz(fi,speci(:,:,i))));
He_sep(:,:,i)=4.*sqrt([He_1;He_2;He_3;He_4;He_5]);
end

% Set up matrix for PCA
matriz=ones(length(He_sep),60)*0;
for i=1:length(He_sep)
matriz(i,:)=reshape(squeeze(He_sep(:,:,i))',1,60);
end
medias=mean(matriz);
std_matriz=std(matriz);

% PCA
limite=90  ;
[coeff,score,latent,tsquared,explained,mu] = pca(matriz);
acum=explained;
for i=2:length(acum)
    acum(i,:)=acum(i,:)+acum(i-1,:);
end
ind_pca_aux=find(acum>limite);
ind_pca=ind_pca_aux(1);

clusters=9;
[idx,C,sumd,D] = kmeans(score(:,1:ind_pca),clusters,'Replicates',20); 
sumD(clusters)=sum(sumd);
meanD(clusters)=mean(sumd);
T_matrix=ones(24,1)*S.T';
T_matrix=T_matrix';

% Recomposição dos espectros
for i=1:clusters
    PCA_matriz(i,:)=C(i,1:ind_pca)*coeff(:,1:ind_pca)';
    PCA_matriz(i,:)=PCA_matriz(i,:)+medias;
    centroid(:,:,i)=(reshape(PCA_matriz(i,:),12,5))';
    [e,r]=find(idx==i);
    ocorr(i,:)=length(e)./length(idx); 
    centroid_original(:,:,i)=mean(S.spec(:,:,e),3);
    Hs_cluster(i,:)=4*sqrt(trapz(S.theta,trapz(S.f,centroid_original(:,:,i))));
    power_centroids(i,:)=density*grav*trapz(S.theta,trapz(S.f,squeeze(centroid_original(:,:,i).*(grav/(4*pi))).*T_matrix));
    m0x=trapz(S.theta,trapz(S.f,centroid_original(:,:,i)));
    m_1x=trapz(S.theta,trapz(S.f,T_matrix.*centroid_original(:,:,i)));
    Tex=m_1x/m0x;
    power_centroids_param(i,:)=(1/(64*pi))*density*grav^2*Hs_cluster(i,:)^2*Tex;
    clear m0x m_1x Tex 
end

% Coloca em ordem crescente de altura
[a,b]=sort(Hs_cluster);
Hs_cluster=Hs_cluster(b,:);
centroid=centroid(:,:,b);
ocorr=ocorr(b,:);
PCA_matriz=PCA_matriz(b,:);
power_centroids=power_centroids(b,:);
centroid_original=centroid_original(:,:,b);
power_centroids_param= power_centroids_param(b,:);


% Plot dos Centroids
for i=1:clusters
h=figure
subplot(1,2,1)
polar_spect_part(5,med_T,centroid(:,:,i));colorbar;colormap(flipud(hot));shading flat
title(['Centroid Cluster #' num2str(i) ' - Hs=' num2str(round(Hs_cluster(i,:),2)) ' - Ocorrência:' num2str(round(ocorr(i,:).*100,2)) '%'])

subplot(1,2,2)
polarcont_naut(S.T,S.theta,centroid_original(:,:,i));colorbar;colormap(flipud(hot));shading flat

set(gcf,'paperposition',[0 0 20 12]);
% print(h,'-dpng',['Centroid_cluster_' num2str(i) '.png']);
end

f_matrix=ones(24,1)*S.f';
f_matrix=f_matrix';
theta_matrix=ones(25,1)*S.theta';
theta_matrix=theta_matrix';

for i=1:length(S.spec)
m0=trapz(S.theta,trapz(S.f,S.spec(:,:,i)));
m_1=trapz(S.theta,trapz(S.f,T_matrix.*S.spec(:,:,i)));
m2=trapz(S.theta,trapz(S.f,f_matrix.^2.*S.spec(:,:,i)));
Tm02(i,:)=sqrt(m0/m2);
Te(i,:)=m_1/m0;
Hm0(i,:)=4.*sqrt(m0);
power_spec(i,:)=density*grav*trapz(S.theta,trapz(S.f,squeeze(S.spec(:,:,i).*(grav/(4*pi))).*T_matrix));
power_param(i,:)=(1/(64*pi))*density*grav^2*Hm0(i,:)^2*Te(i,:);
power_param_periodo_medio(i,:)=(1/(64*pi))*density*grav^2*Hm0(i,:)^2*Tm02(i,:);
f_spec(i,:)=trapz(S.theta,S.spec(:,:,i)');
[~,b]=max(f_spec(i,:));
if b<3
p = polyfit(S.f(b-1:b+1)',f_spec(i,b-1:b+1),2);
x=S.f(b-1):0.001:S.f(b+1); 
elseif b<23
p = polyfit(S.f(b-2:b+2)',f_spec(i,b-2:b+2),2);
x=S.f(b-2):0.001:S.f(b+2);
elseif b<24
p = polyfit(S.f(b-2:b+1)',f_spec(i,b-2:b+1),2);
x=S.f(b-2):0.001:S.f(b+1);
else
p = polyfit(S.f(b-2:b)',f_spec(i,b-2:b),2);
x=S.f(b-2):0.001:S.f(b);
end
[a,b]=max(polyval(p,x));
Tp(i,:)=1/x(b);
A=trapz(S.f,trapz(S.theta,cos(theta_matrix).*S.spec(:,:,i)'));
B=trapz(S.f,trapz(S.theta,sin(theta_matrix).*S.spec(:,:,i)'));
if A<0
    MeanDir(i,:)=rad2deg(atan(B/A)+pi);
else
    MeanDir(i,:)=rad2deg(atan(B/A));
end
end
MeanDir(MeanDir<0)=MeanDir(MeanDir<0)+360;

figure
plot(S.time,Tp,'k',S.time,Te,'r',S.time,Tm02,'b')
legend('Tp','Te','Tm02')
grid on

figure
plot(S.time,power_param,'k',S.time,power_param_periodo_medio,'r')
legend('Te','Tm02')
grid on





%Separa casos Hs-Tp
[classes_hs,classes_te,Hm0_medio,Te_medio,Dir_medio,probabilidade]=selec_casos3(Hm0,Te,MeanDir,power_param,'Condor');

%Compara Energia media gerada com espectros centroides e com casos Hs-Tp

E_real=mean(power_spec)
E_real_param=mean(power_param)
E_centroids=sum(power_centroids.*ocorr)
E_centroids_param=sum(power_centroids_param.*ocorr)

Hm0_medio=reshape(Hm0_medio,25,1);
Te_medio=reshape(Te_medio,25,1);
Dir_medio=reshape(Dir_medio,25,1);

probabilidade=reshape(probabilidade,25,1);
for i=1:length(Te_medio)
power_parametros(i,:)=(1/(64*pi))*density*grav^2*Hm0_medio(i,:)^2*Te_medio(i,:);
end
E_parametros=nansum(power_parametros.*probabilidade*.01)

figure
plot(S.time,power_spec,'k',S.time,power_param,'r');

Results_Clusters.Spectra=S.spec;
Results_Clusters.time=S.time;
Results_Clusters.T=S.T;
Results_Clusters.f=S.f;
Results_Clusters.theta=S.theta;
Results_Clusters.Hm0=Hm0;
Results_Clusters.MeanDir=MeanDir;
Results_Clusters.Tp=Tp;
Results_Clusters.Te=Te;
Results_Clusters.Power_spec=power_spec;
Results_Clusters.Power_param=power_param;
Results_Clusters.Centroids=centroid;
Results_Clusters.Centroids_Originais=centroid_original;
Results_Clusters.Hs_Centroids=Hs_cluster;
Results_Clusters.Power_Centroids_spec=power_centroids;
Results_Clusters.Power_Centroids_param=power_centroids_param;
Results_Clusters.idx=idx;
Results_Clusters.Probabilidade_clusters=ocorr;
Results_Clusters.Casos_Hm0=Hm0_medio;
Results_Clusters.Casos_Te=Te_medio;
Results_Clusters.Casos_MeanDir=Dir_medio;
Results_Clusters.Probabilidades_Casos_Hm0_Te=probabilidade;
save('Condor_Clusters.mat', 'Results_Clusters');

exporta1=[probabilidade Hm0_medio Te_medio Dir_medio];
save('Casos_Hm0-Te.txt','exporta1','-ascii')
[~,ind_hs]=max(Hm0);[~,ind_te]=max(Te);[~,ind_P]=max(power_param);
exporta2=[Hm0(ind_hs) Te(ind_hs) MeanDir(ind_hs);Hm0(ind_P) Te(ind_P) MeanDir(ind_P);Hm0(ind_te) Te(ind_te) MeanDir(ind_te)];
save('Extremos_Hm0-Te.txt','exporta2','-ascii')
for i=1:clusters
    linear_freq(i,:)=trapz(S.theta,centroid_original(:,:,i)');
end
exporta3=[ocorr linear_freq];
save('Casos_Espectros.txt','exporta3','-ascii')
exporta4=[trapz(S.theta,S.spec(:,:,ind_hs)');trapz(S.theta,S.spec(:,:,ind_P)');trapz(S.theta,S.spec(:,:,ind_te)')];
save('Extremos_Espectros.txt','exporta4','-ascii')
