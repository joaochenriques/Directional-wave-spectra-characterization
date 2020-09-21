% Selecao de casos representativos de onda Hs-Tp
% Por Pedro Ribeiro, em 19/03/2018

function [classes_hs,classes_tp,Hs_medio,Tp_medio,Dir_medio,probabilidade]=selec_casos(Hs,Tp,Dir,fluxo,filename)

dhs=((max(Hs))-(min(Hs)))/3;
dtp=((max(Tp)+.5)-(min(Tp)-1.2))/3;

classes_hs=min(Hs):dhs:max(Hs);
classes_hs=[(floor(classes_hs(1)*10))/10 (round(classes_hs(2:3)*10))/10 (ceil(classes_hs(end)*10))/10];
classes_tp=min(Tp)-1.2:dtp:max(Tp)+.5;
classes_tp=[(floor(classes_tp(1)*10))/10 (round(classes_tp(2:3)*10))/10 (ceil(classes_tp(end)*10))/10];
% fluxo=(1/8)*1025*9.81.*(Hs.^2)*(1/2)*(1.56).*Tp;%calcula escalar fluxo de energia p/ cada caso de onda

% primeira distribuiÃ§Ã£o - Media ponderada pelo fluxo de energia em cada classe de Hs-Tp
for i=1:length(classes_hs)-1
    hs_classe=find (Hs>=classes_hs(i) & Hs<classes_hs(i+1));
    for j=1:length(classes_tp)-1
        tp_classe=find(Tp(hs_classe,:)>=classes_tp(j) & Tp(hs_classe,:)<classes_tp(j+1));
        casos=hs_classe(tp_classe);
        fluxo_caso(i,j)=sum(fluxo(casos,:));
        Tp_medio(i,j)=sum(fluxo(casos,:).*Tp(casos,:))./sum(fluxo(casos,:));
        Hs_medio(i,j)=sum(fluxo(casos,:).*Hs(casos,:))./sum(fluxo(casos,:));
        probabilidade(i,j)=100.*length(casos)./length(Hs(~isnan(Hs)));
    end
end

% Faz nova distribuiÃ§Ã£o na classe de maior ocorrência para aumentar
% resolução
[a,b]=max(probabilidade);
[aa,ind2]=max(a);
ind1=b(:,ind2)
if ind1==1
classes_hs=[classes_hs(ind1):(classes_hs(ind1+1)-classes_hs(ind1))/3:classes_hs(ind1+1) classes_hs(ind1+2:end)];    
else
classes_hs=[classes_hs(1:ind1-1) classes_hs(ind1):(classes_hs(ind1+1)-classes_hs(ind1))/3:classes_hs(ind1+1) classes_hs(ind1+2:end)];
end

if ind2==1
classes_tp=[classes_tp(ind2):(classes_tp(ind2+1)-classes_tp(ind2))/3:classes_tp(ind2+1) classes_tp(ind2+2:end)];    
else
classes_tp=[classes_tp(1:ind2-1) classes_tp(ind2):(classes_tp(ind2+1)-classes_tp(ind2))/3:classes_tp(ind2+1) classes_tp(ind2+2:end)];
end

classes_hs=(round(classes_hs*10))/10;
classes_tp=(round(classes_tp*10))/10;

clear probabilidade Tp_medio Hs_medio fluxo_caso

[u,v]=magdir2uv(ones(length(Hs),1),Dir,0);
for i=1:length(classes_hs)-1
    hs_classe=find (Hs>=classes_hs(i) & Hs<classes_hs(i+1));
    for j=1:length(classes_tp)-1
        tp_classe=find(Tp(hs_classe,:)>=classes_tp(j) & Tp(hs_classe,:)<classes_tp(j+1));
        casos=hs_classe(tp_classe);
        fluxo_caso(i,j)=sum(fluxo(casos,:));
        Tp_medio(i,j)=sum(fluxo(casos,:).*Tp(casos,:))./sum(fluxo(casos,:));
        Hs_medio(i,j)=sum(fluxo(casos,:).*Hs(casos,:))./sum(fluxo(casos,:));
        u_medio(i,j)=sum(fluxo(casos,:).*u(casos,:))./sum(fluxo(casos,:));
        v_medio(i,j)=sum(fluxo(casos,:).*v(casos,:))./sum(fluxo(casos,:));
        probabilidade(i,j)=100.*length(casos)./length(Hs(~isnan(Hs)));
    end
end
[mag,Dir_medio]=uv2veldir(u_medio,v_medio);

% Distribuição somente para ilustrar o plot

classes_hs_plot=min(Hs):.3:max(Hs);
classes_tp_plot=min(Tp):.3:max(Tp);

% Media ponderada pelo fluxo de energia em cada classe de Hs-Tp
for i=1:length(classes_hs_plot)-1
    hs_classe_plot=find (Hs>=classes_hs_plot(i) & Hs<classes_hs_plot(i+1));
    for j=1:length(classes_tp_plot)-1
        tp_classe_plot=find(Tp(hs_classe_plot,:)>=classes_tp_plot(j) & Tp(hs_classe_plot,:)<classes_tp_plot(j+1));
        casos=hs_classe_plot(tp_classe_plot);
        fluxo_caso_plot(i,j)=sum(fluxo(casos,:));
        Tp_medio_plot(i,j)=sum(fluxo(casos,:).*Tp(casos,:))./sum(fluxo(casos,:));
        Hs_medio_plot(i,j)=sum(fluxo(casos,:).*Hs(casos,:))./sum(fluxo(casos,:));
        probabilidade_plot(i,j)=100.*length(casos)./length(Hs(~isnan(Hs)));
    end
end


% Plots

ccc=flipud(hot);
ccc=ccc(14:end,:);

probabilidade_plot(probabilidade_plot==0)=nan;

h=figure
pcolor(classes_tp_plot(1:end-1),classes_hs_plot(1:end-1),probabilidade_plot); 
% colormap(flipud(gray));
shading flat; colormap(jet); hold on
hcb=colorbar;
title(hcb,'Occur. (%)')
for i=1:length(classes_hs)-1
plot(Tp_medio(i,:),Hs_medio(i,:),'.k');hold on
end
for i=1:length(classes_hs)
    plot(classes_tp,classes_hs(i)*ones(length(classes_tp),1),'-k');
end
for i=1:length(classes_tp)
    plot(classes_tp(i)*ones(length(classes_hs),1),classes_hs,'-k');
end
hold off
set(gca,'xtick',classes_tp,'xticklabel',classes_tp,'fontsize',6);
set(gca,'ytick',classes_hs,'yticklabel',classes_hs,'fontsize',6);
grid on
xlim([classes_tp(1) classes_tp(end)])
ylim([classes_hs(1) classes_hs(end)])
xlabel('Te (s)');
ylabel('Hm0 (m)');
title('Probability Distribution and Hm0-Te Selected Cases')
set(gcf,'paperposition',[0 0 8 4]);
% print(h,'-dpng','-r300','Casos_Hs_Tp.png')
% close(h)
set(gcf,'paperposition',[0 0 15 10]);
% print(h,'-dpng',[ filename '_Dist_Hs_Te.png']);



% Exporta resultado em TXT
% matriz=[classes_hs' [classes_tp(2:end) ; probabilidade]];
% dlmwrite('probabilidade_Hs_Tp.txt',matriz,'delimiter','\t','precision',3);
% matriz2=[classes_hs' [classes_tp(2:end) ; Hs_medio]];
% dlmwrite('casos_Hs.txt',matriz2,'delimiter','\t','precision',3);
% matriz3=[classes_hs' [classes_tp(2:end) ; Tp_medio]];
% dlmwrite('casos_Tp.txt',matriz3,'delimiter','\t','precision',3);

% save('probabilidades.mat', 'classes_hs','classes_tp','Hs_medio','Tp_medio','probabilidade');
end
