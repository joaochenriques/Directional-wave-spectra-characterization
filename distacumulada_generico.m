function distacumulada_generico(dados, varargin)

%Diagrama de Dsitribuição acumulada com quartis 25%, 50%, 75% e max Dados.
%
% Varargin:

%       - maxdados - valor máximo do eixo x, default é o maior dado da série%
%       - classes: vetor de classes de frequencia. Default é [0 20 40 60 80 100]%
%       - lcb: Legenda do MaxDados. Ex: 'Vel. máx'%
%       - titulo: Título do Gráfico. Ex: 'titulo'
%       - variavel: Legenda do eixo X. Ex: 'Velocidade (nós)'
%       - printf: Salva figura. Igual a zero nao salva. Maior que zero
%         assume o valor de entrada como dpi para salvar em -png.
%       - sufix: nome do arquivo ao salvar figura
%       - Subplot: Plota no axes ativo ou nova figura
%       - parent, por defaultum novo axes é criado, caso contrario parent
%         identifica o axes a ser utilizado

%%
maxdados=max(dados); % valor máximo do eixo x
classes=[0 20 40 60 80 100]; % etor de classes de frequencia
lcb='Vel. máx'; % legenda da barra de cores
titulo=[]; % titulo do gráfico
variavel=[]; % Legenda do eixo X
printf = 0; %salva figura
sufix = ''; % nome ao salvar figura
Subplot = 0; % 0 = nova figura / 1 = plota no axes existente
parent=0; % 0 = novo axes / outro valor identifica o axes existente

%varargin
vin=varargin;
if isempty(vin)==0
    for i=1:length(vin)
        if isequal(vin{i},'maxdados')
            maxdados=vin{i+1};
        elseif isequal(vin{i},'classes')
            classes=vin{i+1};
        elseif isequal(vin{i},'lcb')
            lcb=vin{i+1};
        elseif isequal(vin{i},'titulo')
            titulo=vin{i+1};
        elseif isequal(vin{i},'variavel')
            variavel=vin{i+1};
        elseif isequal(vin{i},'printf')
            printf=vin{i+1};
        elseif isequal(vin{i},'sufix')
            sufix=vin{i+1};
        elseif isequal(vin{i},'Subplot')
            Subplot = vin{i+1};
        elseif isequal(vin{i},'parent')
            parent=vin{i+1};
        end
    end
end

frequencia=linspace(0, max(dados),100);
for i=0:length(frequencia)-1
    if i~0;
        frequencia(2,i)=sum(dados>frequencia(1,i) & dados<=frequencia(1,i+1));
        interv_Hs(i)=(frequencia(1,i)+frequencia(1,i+1))/2;
    end
end
frequencia=(frequencia(2,:)/(sum(frequencia(2,:))))*100;
% bar(interv_Hs(1,:),frequencia(1,1:end-1),1)
frequencia=frequencia';

fr_acum=linspace(0,max(dados),100);
fr_acum=fr_acum';
fr_acum(:,2)=cumsum(frequencia);

if printf > 0
    figure('visible','off');
else
    if Subplot == 1
    else
        figure('visible','on');
    end
end

if parent
  axp=parent;
else
  axp=axes;
end

plot(fr_acum(:,1),fr_acum(:,2), 'r', 'linewidth',2)
axis([0 ceil(max(dados)) 0 100])
grid on
hold on

D=ceil(max(dados))/100;
% encontra o quartil 25
linhah_i25(1,1:101)=0:D:ceil(max(dados));
linhah_i25(2,1:101)=25;
[int25,aux1]=curveintersect(linhah_i25(1,:),linhah_i25(2,:),fr_acum(:,1),fr_acum(:,2));
linhah_i25=[0 int25];
linhah_i25(2,:)=25;
plot(linhah_i25(1,:),linhah_i25(2,:), '-k', 'linewidth',2);
hold on
linhav_i25(2,:)=[int25 int25];
linhav_i25(1,:)=[0 25];
plot(linhav_i25(2,:),linhav_i25(1,:), '-k', 'linewidth',2);

% quartil 50
linhah_i50(1,1:101)=0:D:ceil(max(dados));
linhah_i50(2,1:101)=50;
[int50,aux1]=curveintersect(linhah_i50(1,:),linhah_i50(2,:),fr_acum(:,1),fr_acum(:,2));
linhah_i50=[0 int50];
linhah_i50(2,:)=50;
plot(linhah_i50(1,:),linhah_i50(2,:), '-k', 'linewidth',2);
hold on
linhav_i50(2,:)=[int50 int50];
linhav_i50(1,:)=[0 50];
plot(linhav_i50(2,:),linhav_i50(1,:), '-k', 'linewidth',2);

% quartil 75
linhah_i75(1,1:101)=0:D:ceil(max(dados));
linhah_i75(2,1:101)=75;
[int75,aux1]=curveintersect(linhah_i75(1,:),linhah_i75(2,:),fr_acum(:,1),fr_acum(:,2));
linhah_i75=[0 int75];
linhah_i75(2,:)=75;
plot(linhah_i75(1,:),linhah_i75(2,:), '-k', 'linewidth',2);
hold on
linhav_i75(2,:)=[int75 int75];
linhav_i75(1,:)=[0 75];
plot(linhav_i75(2,:),linhav_i75(1,:), '-k', 'linewidth',2);

plot(fr_acum(:,1),fr_acum(:,2), 'r', 'linewidth',2)
axis([0 ceil(max(dados)) 0 100])
grid on
hold on
set(gca,'fontname','arial','fontsize',8)
xlabel(variavel ,'fontname','arial','fontsize',8)
xlim([0 maxdados]);
ylabel('Cumulative Probability (%)','fontname','arial','fontsize',8)
text(.6*maxdados,3,['Q. 25 % = '  sprintf('%2.2f',int25)],'fontname','arial','fontsize',8)
text(.6*maxdados,8,['Q. 50 % = '  sprintf('%2.2f',int50)],'fontname','arial','fontsize',8)
text(.6*maxdados,13,['Q. 75 % = '  sprintf('%2.2f',int75)],'fontname','arial','fontsize',8)
text(.6*maxdados,18,[lcb ' = '  sprintf('%2.2f',(max(dados)))],'fontname','arial','fontsize',8)
title(titulo ,'fontsize',10,'fontname','arial')

if printf > 0
    set(gcf,'papersize',[8*1.3 8*1.1],'paperorientation','portrait','PaperPosition',[1 1 8*1.3 8*1.1])
    set(gcf,'inverthardcopy','off'); set(gcf,'color',[1 1 1]);
    eval(['print(gcf,''-zbuffer'',''-dpng'',','''-r',num2str(printf),''',''',cd,'\diag_direc',sufix,'.png',''');']);
    close all;
end