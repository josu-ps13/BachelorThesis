%% BACHELOR'S THESIS - Supplementary script
%
% Thesis title: A compartmental model to investigate intracranial pulsatility
% Author: Josue Perez Sabater (josu.ps13@gmail.com)
% Supervisor: Wilfried Coenen (wcoenen@ing.uc3m.es)
%
% Choose a region R and vessel V to visualize the change in velocity along
% time, for each series. Figure shows the velocity field in the chosen ROI
% both with aliasing (for all series) and without aliasing (after
% correcting it). Run this code after running getQ.m

venc=[10 30 80];
ser=1:3;
r=3;
v=1;
auto=1;

lim{1,1}=[106 143;116 157];
lim{2,1}=[ 67  81;103 117];
lim{2,3}=[105 115;122 132];
lim{2,4}=[139 146;131 138];
lim{3,1}=[ 79  90;102 113];
lim{3,2}=[158 168;103 112];
lim{3,3}=[ 96 105;116 126];
lim{3,4}=[144 154;117 126];

col{1,1}=[ -3  3];
col{2,1}=[-25 10];
col{2,3}=[-10 10];
col{2,4}=[ -8 10];
col{3,1}=[  0 40];
col{3,2}=[-50  2];
col{3,3}=[-30  1];
col{3,4}=[-25 .5];

tim{1,1}=1:30;
tim{2,1}=1;
tim{2,3}=1:30;
tim{2,4}=1:30;
tim{3,1}=1:30;
tim{3,2}=1:30;
tim{3,3}=1:30;
tim{3,4}=1:30;

fg=figure(1);
for t=tim{r,v}
    for s=ser
        U_temp=U{r,v}{s}{t};
        U_temp(U_temp==0)=venc(s);
        subplot(2,2,s)
        imshow(U_temp)
        title(SER{s},'fontsize',15)
        clim(venc(s)*[-1 1]);colorbar
        xlim(lim{r,v}(1,:))
        ylim(lim{r,v}(2,:))
    end

    U_temp=U_corrt{r,v}{t};
    U_temp(U_temp==0)=col{r,v}(2);
    subplot(2,2,4)
    imshow(U_temp)
    title('Combined','fontsize',15)

    if auto,clim([min(U_temp(:)) max(U_temp(:))]);colorbar
    else,clim(col{r,v});colorbar,end
    xlim(lim{r,v}(1,:))
    ylim(lim{r,v}(2,:))
    drawnow
end

linkaxes(findall(fg,'type','axes'))
set(findobj('type','fig'),'color','w')
