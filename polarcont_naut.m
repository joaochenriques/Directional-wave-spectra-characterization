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

[a,b] = size(z);

if a ~= length(r)
    error('r is not the same length as the first dimension of z')
end

if b ~= length(theta)
    error('theta is not the same length as the second dimension of z')
end

%Convert to nautical coordinates
thetat=deg2rad(cart2naut(rad2deg(theta)));

x = zeros(a,b);
y = zeros(a,b);

for j = 1:a
    for k = 1:b
        x(j,k) = r(j)*cos(thetat(k));
        y(j,k) = r(j)*sin(thetat(k));
    end
end

% Add column to data to close the circle
x=[x x(:,1)];
y=[y y(:,1)];
z=[z z(:,1)];


pcolor(x,y,z);
axis off
axis image
hold on
plot ([min(min(x)) max(max(x))],[0 0],':k'); hold on
plot ([0 0],[min(min(y)) max(max(y))],':k'); hold on
% text(min(min(x))*1.2,0,'270','fontsize',8);
% text(max(max(x))*1.08,0,'90','fontsize',8);
% text(-1.5,min(min(y))*1.1,'180','fontsize',8);
% text(-1.5,max(max(y))*1.1,'360','fontsize',8);

pos=.25*pi;
th = 0:pi/50:2*pi;
plot(4*cos(th), 4*sin(th),':k');hold on;  text (4.4*cos(pos),-4.4*sin(pos),'4s','fontsize',8);
plot(8*cos(th), 8*sin(th),':k');hold on;  text (8.4*cos(pos),-8.4*sin(pos),'8s','fontsize',8);
plot(12*cos(th), 12*sin(th),':k');hold on;  text (12.4*cos(pos),-12.4*sin(pos),'12s','fontsize',8);
plot(16*cos(th), 16*sin(th),':k');hold on;  text (16.4*cos(pos),-16.4*sin(pos),'16s','fontsize',8);
plot(20*cos(th), 20*sin(th),':k');hold on;  text (20.4*cos(pos),-20.4*sin(pos),'20s','fontsize',8);
end