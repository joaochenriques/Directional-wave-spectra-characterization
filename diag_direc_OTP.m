function diag_direc(mag,dir,varargin)
% diag_direc: DIAGRAMA DIRECIONAL DE FREQUENCIA RELATIVA
%
%   sintax: diag_direc(mag,dir,varargin)
%
%   IDEAL PARA SÉRIES HISTÓRICAS GRANDES, ACIMA DE 30.000 DADOS
%
%   - mag: vetor de magnitudes
%   - dir: vetor de direções (náutica)
%
% Varargin:
%       - nclass: Quantidade de classes de magnitude onde são calculas as
%               frequencias. O valor default é 80, fica mais apresentável.%
%       - ps: Discretização das direções em graus. O valor default é 1, fica
%               mais apresentável. SOMENTE UTILIZAR VALORES INTEIROS DE
%               DIREÇÃO.
%       - eixo: vetor dos círculos de magnitude. Ex: [5 10 15 20]%
%       - dirl: Quadrante CARTESIANO em que será plotada as legendas do
%       - un: unidades dos dados. Ex: 'm/s'%
%               circulos de magnitude. Default é 2%
%       - lcb: Legenda da barra de cores. Ex: 'Frequencia'%
%       - bt: Tamanho do círculo central (0,0). Default é 10%
%       - el: Espessura da linha dos círculos de magnitude. Default é 1%
%       - titulo: Título do Gráfico. Ex: 'titulo'
%       - printf: Salva figura. Igual a zero nao salva. Maior que zero
%       assume o valor de entrada como dpi para salvar em -png.
%       - sufix: nome do arquivo ao salvar figura
%       - quiet (ProgressBar): não mostra a waitbar
%       - Subplot: Plota no axes ativo ou nova figura
%       -parent, por defaultum novo axes é criado, caso contrario parent
%       identifica o axes a ser utilizado

%%
nclass=80; % numero de classes da magnitude
ps=1; % discretização da direção
eixo=[]; % circulos de magnitude
un='m/s'; % unidade dos dados
dirl=2; % quadrante cartesiano da legenda
lcb='Frequência de Ocorrência (%)'; % legenda da barra de cores
bt=10; % tamanho do circulo central 
el=1; % espessura da linha dos circulos de magnitude
titulo=[]; % titulo do gráfico
printf = 0; %salva figura
sufix = ''; % nome ao salvar figura
ProgressBar = '';
Subplot = 0; % 0 = nova figura / 1 = plota no axes existente
parent=0; % 0 = novo axes / outro valor identifica o axes existente

% offset da legenda dos circulos de %
if nanmax(mag)>16
    diff=1;
elseif nanmax(mag) > 3
    diff=0.5;
elseif nanmax(mag) > 1 &  nanmax(mag) <=3
    diff=0.25;
elseif nanmax(mag) <= 1
    diff=0.1;
elseif nanmax(mag) <= 0.5
    diff=0.05;
end

%varargin
vin=varargin;
if isempty(vin)==0
    for i=1:length(vin)
        if isequal(vin{i},'nclass')
            nclass=vin{i+1};
        elseif isequal(vin{i},'ps')
            ps=vin{i+1};
        elseif isequal(vin{i},'eixo')
            eixo=vin{i+1};
        elseif isequal(vin{i},'un')
            un=vin{i+1};
        elseif isequal(vin{i},'dirl')
            dirl=vin{i+1};
        elseif isequal(vin{i},'lcb')
            lcb=vin{i+1};
        elseif isequal(vin{i},'bt')
            bt=vin{i+1};
        elseif isequal(vin{i},'el')
            el=vin{i+1};
        elseif isequal(vin{i},'titulo')
            titulo=vin{i+1};
        elseif isequal(vin{i},'printf')
            printf=vin{i+1};
        elseif isequal(vin{i},'sufix')
            sufix=vin{i+1};
        elseif isequal(vin{i},'quiet')
            ProgressBar='quiet';
        elseif isequal(vin{i},'Subplot')
            Subplot = vin{i+1};
        elseif isequal(vin{i},'parent')
            parent=vin{i+1};
        end
    end
end

% Definições
if isempty(eixo)==1
    if nanmax(mag) > 10 % vento
        aux01 = nanmax(mag); aux01 = (aux01+(aux01*.2))/4;
        eixo = round([aux01 aux01*2 aux01*3 aux01*4]);
    else % onda
        if nanmax(mag) > 3
            eixo = 1 : 1 : ceil(nanmax(mag));
        elseif nanmax(mag) > 1 & nanmax(mag) <= 3
            eixo = 0 : 0.75 : ceil(nanmax(mag));
        elseif nanmax(mag) <= 1
            eixo = 0 : 0.25 : ceil(nanmax(mag));
        end
%         eixo=[ceil(nanmin(mag)+nanmin(mag)*0.7) ceil(nanmedian(mag)+nanmedian(mag)*0.7) ceil(nanmax(mag))];
    end
end

if isempty(titulo)==1
    %titulo='Diagrama Direcional de Freqûencia Relativa';
    titulo = '';
end

% calculo de parametros
% interv=(max(mag(:,1)) - min(mag(:,1)))/nclass;
classes=linspace(min(mag(:,1)),max(mag(:,1)),nclass);%[min(mag(:,1)) min(mag(:,1))+interv : interv : max(mag(:,1))];    
interv=classes(2)-classes(1);
classes=[classes classes(end)+interv];

D=max(eixo);
quad=((dirl/4)/0.5)*pi - 67*(pi/180);
for g=1:length(eixo)
    textL{g,1}=[num2str(eixo(1,g)),' ',un];
end
piC=(180/ps)-1;

% classes de referencia
dir0=0:ps:360; dir0=dir0';

if printf > 0
    figure('visible','off');
else
    if Subplot == 1
    else
        figure('visible','on');
    end
end

% Cálculo das frequências nas classes
cont=zeros(length(dir0)-1,nclass);
if ~strcmp(ProgressBar,'quiet')
    O = waitbar(0,'Espere... calculando a frequência nas classes');
end
g0 = 1;
for h=1:nclass
    for g=2:length(dir0)
        c=find(mag > classes(1,h) & mag <= classes(1,h+1));
        d=find(dir(c,1) > dir0(g-1,1) & dir(c,1) <= dir0(g,1));
        cont(g-1,h)=length(d);
        clear c d
    end
    if ~strcmp(ProgressBar,'quiet')
        if g0 == 1
            waitbar(0,O,'Espere... calculando a frequência nas classes .') ; g0 = g0 + 1;
        elseif g0 == 2
            waitbar(0,O,'Espere... calculando a frequência nas classes ..') ; g0 = g0 + 1;
        elseif g0 == 3
            waitbar(0,O,'Espere... calculando a frequência nas classes ...') ; g0 = g0 + 1;
        elseif g0 == 4
            waitbar(0,O,'Espere... calculando a frequência nas classes .') ; g0 = 1;
        end        
        waitbar((h)/nclass,O); drawnow;
    end
end
if ~strcmp(ProgressBar,'quiet')
    close(O)
end
cont(cont==0)=nan; cont=(cont*100)/length(mag);
%%
% cria grid circular
for h=1:nclass
    r=classes(1,h+1) - (interv/2);
    t=linspace(pi/2,2*pi + pi/2,length(dir0)-1);%pi/2 : pi/piC : 2*pi + pi/2;
    x(:,h)=r*cos(t);
    y(:,h)=r*sin(t);
end

% circulos de magnitude
for g=1:length(eixo)
    r=eixo(1,g);
    t=pi/2 : pi/piC : 2*pi + pi/2;
    xE(:,g)=r*cos(t);
    yE(:,g)=r*sin(t);
end

% legenda dos circulos de magnitude
for g=1:length(eixo)
    r=eixo(1,g) - diff;
    if dirl==1
        xL(:,g)= (r*sin(quad));
        yL(:,g)= (r*cos(quad));
    elseif dirl==2        
        xL(:,g)=(r*cos(quad));
        yL(:,g)=(r*sin(quad));
    elseif dirl==3
        xL(:,g) = (r*sin(quad));
        yL(:,g) = (r*cos(quad));
    elseif dirl==4
        xL(:,g) = - (r*sin(quad));
        yL(:,g) = - (r*cos(quad));
    end
end

% plot
if parent
  axp=parent;
else
  axp=axes;
end

pcolor(-x,y,cont(1:end,:)); shading flat; axis equal; hold on
plot(xE,yE,'k:','linewidth',el)
xlim([-D-(D/10*9) D+(D/10*5)]); ylim([-D-(D/10*8) D+(D/10*2)]);

% eixos N,S,E,W e outros
plot([0 D],[0 0],'k:'); text(D+D*0.05,0,'E','fontsize',8) %eixo leste
plot([0 0],[0 D],'k:'); text(-D*0.015,D+D*0.05,'N','fontsize',8)% eixo norte
plot([-D 0],[0 0],'k:'); text(-D-D*0.07,0,'W','fontsize',8)%eixo oeste
plot([0 0],[-D 0],'k:'); text(-D*0.015,-D-D*0.05,'S','fontsize',8)%eixo sul
% plot([0 -D*cos(135*pi/180)],[0 -D*sin(135*pi/180)],'k-'); text(-D*cos(135*pi/180),-D*sin(135*pi/180)-D*0.05,'SE','fontsize',8)
% plot([0 D*cos(45*pi/180)],[0 D*sin(45*pi/180)],'k-'); text(D*cos(45*pi/180)+D*0.01,D*sin(45*pi/180)+D*0.05,'NE','fontsize',8)
% plot([0 -D*cos(45*pi/180)],[0 -D*sin(45*pi/180)],'k-'); text(-D*cos(45*pi/180)-D*0.08,-D*sin(45*pi/180)-D*0.05,'SW','fontsize',8)
% plot([0 D*cos(135*pi/180)],[0 D*sin(135*pi/180)],'k-'); text(D*cos(135*pi/180)-D*0.08,D*sin(45*pi/180)+D*0.05,'NW','fontsize',8)
% 
% configurações do plot

axpp=axp;
set(gcf,'currentaxes',axpp); box off
plot(0,0,'ko','markersize',bt,'markerfacecolor','w'); axis equal
xlim([-D-(D/10*9) D+(D/10*5)]); ylim([-D-(D/10*8) D+(D/10*2)]);
set(axpp,'color','none')
set(axpp,'ytick',[],'yticklabel',[],'xtick',[],'xticklabel',[])
set(axpp,'xcolor','w','ycolor','w')

load('jet_claro3.mat')
text(xL,yL,textL,'fontsize',8,'backgroundcolor','w','horizontalalignment','center')
set(axp,'ytick',[],'yticklabel',[],'xtick',[],'xticklabel',[])
set(gcf,'currentaxes',axp); box off
ccc=colorbar;set(ccc,'Fontsize',8,'location','south','position',[0.41 0.21 0.31 0.02]);
text(0,-D-(D*0.2),lcb,'fontsize',10,'horizontalalignment','center')
set(axp,'xcolor','w'); set(axp,'ycolor','w')
colormap(axp,jet_claro3)
title(titulo,'fontsize',12)
% linkaxes([axp axpp],'xy')
%
% salva figura
if printf > 0
    set(gcf,'papersize',[5*1.3 5*1.1],'paperorientation','portrait','PaperPosition',[1 1 8*1.3 8*1.1])
    set(gcf,'inverthardcopy','off'); set(gcf,'color',[1 1 1]);
    eval(['print(gcf,''-zbuffer'',''-dpng'',','''-r',num2str(printf),''',''',cd,'\diag_direc',sufix,'.png',''');']);
    close all;
end