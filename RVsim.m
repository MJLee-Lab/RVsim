function [RV_timecourse] = RVsim(untreated_tau,treated_tau,treated_dr,finaltime,dose_curve_tp,ec50)

arguments
    untreated_tau
    treated_tau
    treated_dr
    finaltime
    dose_curve_tp
    ec50 = 4
end
    
RV_timecourse = table();

% Initial parameters
init = 100;
drc = 0;
xfit = linspace(0,finaltime,500)';

% Control (live)
LiveCellsControl = init*(2.^(untreated_tau * finaltime));
DeadCellsControl = init*(2.^(untreated_tau * finaltime))*drc;

% Untreated
tc = untreated_tau;
drd = 0;

% End point
qtable = table();
qtable.Treatment(1) = {'Untreated'};
qtable.Death(1) = {'-'};

LiveCellsDrug = init*(2.^(finaltime.*tc)) - init*(2.^(finaltime.*tc))*drd;
DeadCellsDrug = init*(2.^(finaltime.*tc))*drd;
qtable.GrowthRate(1) = tc;
qtable.DeathRate(1) = drd;

% Fit
x = [0,;finaltime];
y = [init;LiveCellsDrug];
qfit = fit(x,y,'exp1');
yfit = qfit(xfit);
qtable.PopulationSize{1} = array2table([xfit,yfit],'VariableNames',[{'Time'},{'LiveCells'}]); %qtable.PopulationSize{1} = array2table([xfit,yfit],'VariableNames',[{'Time'},{'LiveCells'}]);

RV_timecourse = [RV_timecourse; qtable];

clear LiveCellsDrug DeadCellsDrug qtable tc drd x y qfit yfit


% Treated, - death
tc = treated_tau; 
drd = 0;

% End point
qtable = table();
qtable.Treatment(1) = {'Treated'};
qtable.Death(1) = {'-'};

LiveCellsDrug = init*(2.^(finaltime.*tc)) - init*(2.^(finaltime.*tc))*drd;
DeadCellsDrug = init*(2.^(finaltime.*tc))*drd;
qtable.GrowthRate(1) = tc;
qtable.DeathRate(1) = drd;

% Fit
x = [0,;finaltime];
y = [init;LiveCellsDrug];
qfit = fit(x,y,'exp1');
yfit = qfit(xfit);
qtable.PopulationSize{1} = array2table([xfit,yfit],'VariableNames',[{'Time'},{'LiveCells'}]);

RV_timecourse = [RV_timecourse; qtable];

clear LiveCellsDrug DeadCellsDrug qtable tc drd x y qfit yfit


% Treated, + death
tc = treated_tau;
drd = treated_dr;

% End point
qtable = table();
qtable.Treatment(1) = {'Treated'};
qtable.Death(1) = {'+'};

LiveCellsDrug = init*(2.^(finaltime.*tc)) - init*(2.^(finaltime.*tc))*drd;
DeadCellsDrug = init*(2.^(finaltime.*tc))*drd;
qtable.GrowthRate(1) = tc;
qtable.DeathRate(1) = drd;

% Fit
x = [0,;finaltime];
y = [init;LiveCellsDrug];
qfit = fit(x,y,'exp1');
yfit = qfit(xfit);
qtable.PopulationSize{1} = array2table([xfit,yfit],'VariableNames',[{'Time'},{'LiveCells'}]);

RV_timecourse = [RV_timecourse; qtable];

clear LiveCellsDrug DeadCellsDrug qtable tc drd x y qfit xfit yfit

clear init drc DeadCellsControl LiveCellsControl



% Simulation - for different death rates
mtable = table();

init = 100;
drc = 0;
xfit = linspace(0,finaltime,500)';
qtreated_dr = linspace(0.05,.95,19);

% Control (live)
LiveCellsControl = init*(2.^(untreated_tau * finaltime));
DeadCellsControl = init*(2.^(untreated_tau * finaltime))*drc;

% Untreated
tc = untreated_tau;
drd = 0;

% End point
qtable = table();
qtable.Treatment(1) = {'Untreated'};
qtable.Death(1) = {'-'};

LiveCellsDrug = init*(2.^(finaltime.*tc)) - init*(2.^(finaltime.*tc))*drd;
DeadCellsDrug = init*(2.^(finaltime.*tc))*drd;
qtable.GrowthRate(1) = tc;
qtable.DeathRate(1) = drd;

% Fit
x = [0,;finaltime];
y = [init;LiveCellsDrug];
qfit = fit(x,y,'exp1');
yfit = qfit(xfit);
qtable.PopulationSize{1} = array2table([xfit,yfit],'VariableNames',[{'Time'},{'LiveCells'}]);

mtable = [mtable; qtable];

clear LiveCellsDrug DeadCellsDrug qtable tc drd x y qfit yfit

for i = 1:length(qtreated_dr)

    % Treated, - death
    tc = treated_tau;
    drd = 0;

    % End point
    qtable = table();
    qtable.Treatment(1) = {'Treated'};
    qtable.Death(1) = {'-'};

    LiveCellsDrug = init*(2.^(finaltime.*tc)) - init*(2.^(finaltime.*tc))*drd;
    DeadCellsDrug = init*(2.^(finaltime.*tc))*drd;
	qtable.GrowthRate(1) = tc;
    qtable.DeathRate(1) = drd;

    % Fit
    x = [0,;finaltime];
    y = [init;LiveCellsDrug];
    qfit = fit(x,y,'exp1');
    yfit = qfit(xfit);
    qtable.PopulationSize{1} = array2table([xfit,yfit],'VariableNames',[{'Time'},{'LiveCells'}]);

    mtable = [mtable; qtable];

    clear LiveCellsDrug DeadCellsDrug qtable tc drd x y qfit yfit


    % Treated, + death
    tc = treated_tau; 
    drd = qtreated_dr(i);

    % End point
    qtable = table();
    qtable.Treatment(1) = {'Treated'};
    qtable.Death(1) = {'+'};

    LiveCellsDrug = init*(2.^(finaltime.*tc)) - init*(2.^(finaltime.*tc))*drd;
    DeadCellsDrug = init*(2.^(finaltime.*tc))*drd;
    qtable.GrowthRate(1) = tc;
    qtable.DeathRate(1) = drd;

    % Fit
    x = [0,;finaltime];
    y = [init;LiveCellsDrug];
    qfit = fit(x,y,'exp1');
    yfit = qfit(xfit);
    qtable.PopulationSize{1} = array2table([xfit,yfit],'VariableNames',[{'Time'},{'LiveCells'}]);

    mtable = [mtable; qtable];
    
    clear LiveCellsDrug DeadCellsDrug qtable tc drd x y qfit yfit
end
clear i xfit

clear init drc qtreated_dr DeadCellsControl LiveCellsControl


% Plot sensitivity (delta)
uni_death = unique(mtable.DeathRate);

qcolor = colorGradient([217, 220, 222]./255, [161 184 195]./255, 11);
qcolor(11:20,:) = colorGradient([161 184 195]./255, [56, 90, 112]./255, 10);

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',3);

% Range of death rates
for i = 2:length(uni_death)
    
    controltable = mtable((mtable.DeathRate==0 &(strcmp(mtable.Treatment,'Untreated'))),:);
    atable = mtable((mtable.DeathRate==0 &(strcmp(mtable.Treatment,'Treated'))),:);
    btable = mtable((mtable.DeathRate==uni_death(i) &(strcmp(mtable.Treatment,'Treated'))),:);
 
    control = controltable.PopulationSize{1};
    q = atable.PopulationSize{1};  
    rv = (q.LiveCells./control.LiveCells);

    control = controltable.PopulationSize{1};
    q = btable.PopulationSize{1};  
    rv = rv - ((q.LiveCells./control.LiveCells));

    plot(q.Time,rv,'LineStyle','-','Color',qcolor(i,:),'LineWidth',3.5,'HandleVisibility','off')
    hold on
    
    clear controltable atable btable control q rv
end

% Input death rate
controltable = RV_timecourse((RV_timecourse.DeathRate==0 &(strcmp(RV_timecourse.Treatment,'Untreated'))),:);
atable = RV_timecourse((strcmp(RV_timecourse.Death,'-') & (strcmp(RV_timecourse.Treatment,'Treated'))),:);
btable = RV_timecourse((strcmp(RV_timecourse.Death,'+') & (strcmp(RV_timecourse.Treatment,'Treated'))),:);

control = controltable.PopulationSize{1};
q = atable.PopulationSize{1};  
rv = (q.LiveCells./control.LiveCells);
control = controltable.PopulationSize{1};
q = btable.PopulationSize{1};  
rv = rv - ((q.LiveCells./control.LiveCells));

plot(q.Time,rv,'LineStyle','-','Color',[214, 4, 31]./255,'LineWidth',3.5)
hold on
plot([0 168],[0 0],'LineStyle','-','Color','k','LineWidth',2.5,'HandleVisibility','off')
   
xlabel('doublings') 
ylabel('Relative sensitivity to cell death')
xticks([0,24,48,72,96,120,144,168])
xticklabels([0,1,2,3,4,5,6,7])
xlim([0 168])
% ylim([0 .2])
pbaspect([1 1 1])
legend(strcat(string(treated_dr),{' '},'death rate'))
clear q x y figure1 axes1 rv i qcolor death arrest uni_death mtable controltable atable btable control q rv



% Fit dose curves
for i = 1:size(RV_timecourse,1)

    control = RV_timecourse.PopulationSize{1};
    q = RV_timecourse.PopulationSize{i};
    
    rv = q.LiveCells./control.LiveCells;
    
    
    for ii = 1:length(rv)
        xfit = linspace(1,8,100);
        yfit = rv(ii) + (1 - rv(ii)) ./ (1+(10.^((xfit-ec50).*1)));

        time(ii,1) = q.Time(ii);
        qfit{ii,1} = array2table([xfit',yfit'],'VariableNames',[{'Dose'},{'RV'}]);
 
        clear xfit yfit
    end
	RVFit = table();
    RVFit.Time = time;
    RVFit.Fit = qfit;
    
    RV_timecourse.RV_fits{i} = RVFit;
    
    clear control q rv RVFit ii time qfit
end
clear i


% Plot dose curve
c = [175,175,175; 214, 4, 31; 214, 4, 31]./255;
c2 = [200,200,200; 250, 160, 173; 250, 160, 173]./255;
l = [{'-'}; {'--'}; {'-'}];

time = dose_curve_tp;

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',20,'LineWidth',3);

for i = 1:size(RV_timecourse,1)
   
    q = RV_timecourse.RV_fits{i};
    q = q.Fit{(find(abs(q.Time - time) == min(abs(q.Time - time))))};
    
    if strcmp(RV_timecourse.Death{i},'+')
    	plot(q.Dose,q.RV,'LineStyle',l{i},'Color',c2(i,:),'LineWidth',8,'HandleVisibility','off')  
        hold on
    end
    plot(q.Dose,q.RV,'LineStyle',l{i},'Color',c(i,:),'LineWidth',4)
    hold on

    clear q
end
clear i 

xlabel('Dose')
xticks([1:8])
ylabel('RV')
xlim([1 8])
ylim([0 1.05])
pbaspect([1 1 1])

legend('Untreated','Treated (-Death)','Treated (+Death)','Location','SouthWest')
title(strcat(string(time),'hr'))

clear figure1 axes1 time c c2 l 


end

        
%% FUNCTIONS
function [grad,im]=colorGradient(c1,c2,depth)
% COLORGRADIENT allows you to generate a gradient between 2 given colors,
% that can be used as colormap in your figures.
%
% USAGE:
%
% [grad,im]=getGradient(c1,c2,depth)
%
% INPUT:
% - c1: color vector given as Intensity or RGB color. Initial value.
% - c2: same as c1. This is the final value of the gradient.
% - depth: number of colors or elements of the gradient.
%
% OUTPUT:
% - grad: a matrix of depth*3 elements containing colormap (or gradient).
% - im: a depth*20*3 RGB image that can be used to display the result.
%
% EXAMPLES:
% grad=colorGradient([1 0 0],[0.5 0.8 1],128);
% surf(peaks)
% colormap(grad);
%
% --------------------
% [grad,im]=colorGradient([1 0 0],[0.5 0.8 1],128);
% image(im); %display an image with the color gradient.

% Copyright 2011. Jose Maria Garcia-Valdecasas Bernal
% v:1.0 22 May 2011. Initial release.

%Check input arguments.
%input arguments must be 2 or 3.
error(nargchk(2, 3, nargin));

%If c1 or c2 is not a valid RGB vector return an error.
if numel(c1)~=3
    error('color c1 is not a valir RGB vector');
end
if numel(c2)~=3
    error('color c2 is not a valir RGB vector');
end

if max(c1)>1&&max(c1)<=255
    %warn if RGB values are given instead of Intensity values. Convert and
    %keep procesing.
    warning('color c1 is not given as intensity values. Trying to convert');
    c1=c1./255;
elseif max(c1)>255||min(c1)<0
    error('C1 RGB values are not valid.')
end

if max(c2)>1&&max(c2)<=255
    %warn if RGB values are given instead of Intensity values. Convert and
    %keep procesing.
    warning('color c2 is not given as intensity values. Trying to convert');
    c2=c2./255;
elseif max(c2)>255||min(c2)<0
    error('C2 RGB values are not valid.')
end
%default depth is 64 colors. Just in case we did not define that argument.
if nargin < 3
    depth=64;
end

%determine increment step for each color channel.
dr=(c2(1)-c1(1))/(depth-1);
dg=(c2(2)-c1(2))/(depth-1);
db=(c2(3)-c1(3))/(depth-1);

%initialize gradient matrix.
grad=zeros(depth,3);
%initialize matrix for each color. Needed for the image. Size 20*depth.
r=zeros(20,depth);
g=zeros(20,depth);
b=zeros(20,depth);
%for each color step, increase/reduce the value of Intensity data.
for j=1:depth
    grad(j,1)=c1(1)+dr*(j-1);
    grad(j,2)=c1(2)+dg*(j-1);
    grad(j,3)=c1(3)+db*(j-1);
    r(:,j)=grad(j,1);
    g(:,j)=grad(j,2);
    b(:,j)=grad(j,3);
end

%merge R G B matrix and obtain our image.
im=cat(3,r,g,b);
end