% POLARCONT Polar contour plot
%
% Richard Rieber
% rrieber@gmail.com
% April 4, 2007
% Updated June 15, 2007
% 
% function [C,h] = polarcont(r,theta,z,N,s)
%
% Purpose: This function creates polar contour plots on the current active
%          figure
% 
% Inputs:  o r     - Radius vector of length m
%          o theta - Angle vector in radians of length n
%          o z     - Magnitude at the points specified in r and theta of
%                    size m x n
%          o N     - The number of contours to plot [OPTIONAL]
%          o s     - Linespec as described in PLOT [OPTIONAL]
%
% Outputs: o C     - returns contour matrix C as described in CONTOURC
%          o h     - Column vector H of handles to LINE or PATCH objects,
%                    one handle per line.  
%
% OTHER NOTES:
% - Both C and h can be used as inputs to CLABEL
% - Colors are defined in colormap
% - Treat this function as a standard contour plot

function [C,h] = polarcont_naut(r,theta,z)
z=flipud(z);
rr=1:r+1;
for i=2:2:2*length(rr)
    r_aux(i:i+1)=[rr(i/2) rr(i/2)];
end 
rr=r_aux(2:end);

for i=2:2:2*length(theta)
    theta_aux(i:i+1)=[theta(i/2) theta(i/2)];
end 
thetat=theta_aux(2:end);

 for i=2:2:2*r
    z_aux(i:i+1,:)=[z(i/2,:);z(i/2,:)];
 end
z=z_aux(2:end,:); 
clear z_aux
%  for i=2:2:2*length(theta)
%     z_aux(:,i:i+1)=[z(:,i/2) z(:,i/2)];
%  end
% %  z=z_aux(:,2:end);
%  z=z_aux;
% [a,b] = size(z);
% 
% if a ~= length(r)-2
%     error('r is not the same length as the first dimension of z')
% end
% 
% if b ~= length(theta)
%     error('theta is not the same length as the second dimension of z')
% end

%Convert to nautical coordinates
thetat=deg2rad(cart2naut(rad2deg(thetat)));

x = zeros(length(rr),length(thetat));
y = zeros(length(rr),length(thetat));

for j = 1:length(rr)
    for k = 1:length(thetat)
        x(j,k) = rr(j)*cos(thetat(k));
        y(j,k) = rr(j)*sin(thetat(k));
    end
end

% Add column to data to close the circle
x=[x x(:,1)];
y=[y y(:,1)];
z=[z z(:,1)];


for i=1:r
    for j=1:length(theta)
pcolor(x(2*i:2*i+1,2*j:2*j+1),y(2*i:2*i+1,2*j:2*j+1),ones(2,2)*mean(z(2*i,j:j+1)));hold on
    end
end
hold off
axis off
caxis([0 max(max(z))])
axis image
hold on
plot ([min(min(x)) max(max(x))],[0 0],':k'); hold on
plot ([0 0],[min(min(y)) max(max(y))],':k'); hold on
% text(min(min(x))*1.2,0,'270','fontsize',8);
% text(max(max(x))*1.08,0,'90','fontsize',8);
% text(-.2,min(min(y))*1.1,'180','fontsize',8);
% text(-.2,max(max(y))*1.1,'360','fontsize',8);

% text (1.2,.2,'Wind','fontsize',8);
% hold on;  text (2.2,.2,'Inter','fontsize',8);
% hold on;  text (3.2,.2,'Swell','fontsize',8);
% pos=.25*pi;
% th = 0:pi/50:2*pi;
% plot(4*cos(th), 4*sin(th),':k');hold on;  text (4.4*cos(pos),-4.4*sin(pos),'4s','fontsize',8);
% plot(8*cos(th), 8*sin(th),':k');hold on;  text (8.4*cos(pos),-8.4*sin(pos),'8s','fontsize',8);
% plot(12*cos(th), 12*sin(th),':k');hold on;  text (12.4*cos(pos),-12.4*sin(pos),'12s','fontsize',8);
% plot(16*cos(th), 16*sin(th),':k');hold on;  text (16.4*cos(pos),-16.4*sin(pos),'16s','fontsize',8);
% plot(20*cos(th), 20*sin(th),':k');hold on;  text (20.4*cos(pos),-20.4*sin(pos),'20s','fontsize',8);

end