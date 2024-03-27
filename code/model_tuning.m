%% BACHELOR'S THESIS - Script 3/3: Model tuning
%
% Thesis title: A compartmental model to investigate intracranial pulsatility
% Author: Josue Perez Sabater (josu.ps13@gmail.com)
% Supervisor: Wilfried Coenen (wcoenen@ing.uc3m.es)
% 
% This script runs the model simulation using several values of a chosen
% parameter PARAM in order to study its effect on a specific signal (flow
% through SAS, Q_sas). The model used is model C from Script 2/3. Both the
% experimental (original) and simulated signals are plotted together, for
% all values of PARAM.
%
% If a different signal than Q_sas is to be studied, modify the code within
% the 'for' loop below.

close all;clear;clc
load ../OUT/Qn.mat  % flow coefficients (ml/s)
load ../OUT/Pn.mat  % carotid pressure coefficients (mmHg)
Qn_expe.art=Qn{3}'; % experimental signals coefficients
Qn_expe.ven=Qn{2}';
Qn_expe.sas=Qn{1}';

%% Define variables and constants
[Pn,Vn,Qn,Vo,R,C,Pjug,Pout]=Define_VarAndCons(Pn,Qn);

% Data to reconstruct signals
T=2;            % total time (s)
t_simu=0:.01:T; % time (s) (simulated signals)
t_expe=0:.03:T; % time (s) (experimental signals)
f=1;            % cardiac frequency (1/s)
w=2*pi*f;       % angular frequency (rad/s)
N=9;            % number of coefficients
n=(0:N-1)';     % coefficient index (0 to 8)
blue=[.1922 .5098 .7412]; % blue color for plot
Q_expe=Reconstruct(Qn_expe,t_expe,w,n);

% Parameter to be study
param={'C' 'art'};
value=[.01 .1 1 10 100];

%% Input to the model
input={'Q','art'};

% Solve matrix singularization if needed
singul.activate=1; % set to 1 in order to impose an additional constrain
if singul.activate==1,singul.input={'Qsas',0};end % average flow through SAS is zero

%% Run models, reconstruct simulated signals and plot results
VAL=numel(value);       % number of values to be tested
SAD=zeros(size(value)); % sum of absolute differences between Q_EXPE and Q_SIMU

figure
for v=1:VAL
    disp("Running model "+v+"/"+VAL)

    % Define the parameter to study
    switch param{1}
        case 'C',C.(param{2})=value(v);
        case 'R',R.(param{2})=value(v);
        case 'Vo',Vo.(param{2})=value(v);
        case 'P'
            if strcmp(param{2},'jug'),Pjug=value(v);
            else,Pout=value(v);end
    end

    % Run model and reconstruct signals in time domain. To study a signal
    % other than Q_sas (the one used here), modify which type of signal
    % (pressure, volume or flow) returns RUNMODEL and change the specific
    % variable "sas" in the plot.
    [~,~,Qn_simu]=RunModel(Pn,Vn,Qn,Vo,R,C,w,N,Pjug,Pout,input,singul);
    Q_simu=Reconstruct(Qn_simu,t_simu,w,n);

    % Plot experimental and simulated signals of Q_sas
    subplot(1,VAL,v),hold on,grid on
    plot(t_expe,Q_expe.sas,'o-',color=blue,linewidth=1.5,markerfacecolor=blue,markersize=3)
    plot(t_simu,Q_simu.sas,'r',linewidth=1.5)
    title(param{1}+"_{"+param{2}+"} = "+value(v))
    if v==1,ylabel("CSF flow (cm^3/s)"),end
    if v==5,legend({'Experimental signal' 'Simulated signal'}),end
    set(gca,fontsize=15,xticklabels=[])

    SAD(v)=sum(abs(Q_simu.sas(1:3:end)-Q_expe.sas));
end

for v=1:VAL,disp("SAD for "+value(v)+": "+SAD(v)),end
set(findobj('type','fig'),'color','w')

%% %%%%%%%%%%%%%%%%%%% FUNCTIONS USED IN THIS SCRIPT %%%%%%%%%%%%%%%%%%%%%%

function[P,V,Q,Vo,R,C,Pjug,Pout]=Define_VarAndCons(Pn,Qn)
% Pressures (mm Hg)
P.car=Pn; % pressure in carotid artery
P.art=[]; % pressure in arteries
P.ven=[]; % pressure in veins
P.cra=[]; % pressure in cranial SAS
P.spi=[]; % pressure in spinal SAS
P.bra=[]; % pressure in brain

% Volumes (ml)
V.art=[]; % volume of arteries
V.ven=[]; % volume of veins
V.cra=[]; % volume of cranial SAS
V.spi=[]; % volume of spinal SAS

% Flows (ml/s)
Q.art=Qn{3}; % flow through arteries
Q.cap=[];    % flow through capillaries
Q.ven=Qn{2}; % flow through veins
Q.sas=Qn{1}; % flow from cranial to spinal SAS

% Compartment volume at rest (ml)
Vo.art=30;   % volume of arteries
Vo.ven=60;   % volume of veins
Vo.cra=47.5; % volume of cranial SAS
Vo.spi=95;   % volume of spinal SAS

% Compartment compliances (ml/mmHg)
C.art=.025; % compliance of arteries
C.cap=.005; % compliance of capillaries
C.ven=.04;  % compliance of veins
C.cra=.725; % compliance of cranial SAS
C.spi=.022; % compliance of spinal SAS
C.bra=40;   % compliance of brain

% Vessel resistances (mmHg s/ml)
R.art=1.4;      % resistance in carotid arteries
R.ven=.4;       % resistance in jugular veins
R.sas=3.75;     % resistance in cranial-spinal passage
R.cap=1056;     % resistance in capillaries

% Other constants
Pjug=3;  % pressure in jugular vein
Pout=10; % pressure out of spinal SAS
end

function[P,V,Q]=RunModel(P,V,Q,Vo,R,C,w,N,Pjug,Pout,input,singul)

Vtot=sum(cell2mat(struct2cell(Vo)))-Vo.spi; % total cranial volume

% Find position of input unknown
PVQ={P V Q};
X=PVQ{input{1}=='PVQ'}; % whether the unknown is P, V or Q
unknowns=[strcat('P',fieldnames(P));strcat('V',fieldnames(V));strcat('Q',fieldnames(Q))]';
i=contains(unknowns,cell2mat(input)); % position of input

% Run model for each mode
for n=0:N-1 % 0 to 8
    inw=1i*n*w;
    inputA=zeros(1,numel(i));inputA(i)=1; % input coefficients
    inputB=X.(input{2})(n+1);

    % Avoid matrix singularization
    if singul.activate && n==0
        j=contains(unknowns,singul.input{1}); % position of additional input
        inputA2=zeros(1,numel(i));inputA2(j)=1;
        inputB2=singul.input{2};
        inputA=[inputA;inputA2];
        inputB=[inputB inputB2];end

    % Define system AX=B
    A=[... % coefficients of the system
        1 -1  0  0 0 0 0 0 0 0 -R.art 0   0    0;...
        0  1 -1  0 0 0 0 0 0 0   0 -R.cap 0    0;...
        0  0 -1  0 0 0 0 0 0 0   0    0 -R.ven 0;...
        0  0  0 -1 1 0 0 0 0 0   0    0   0 -R.sas;...

        0 C.art 0   0   0 -C.art -1  0  0  0 0 0 0 0;...
        0   0 C.ven 0   0 -C.ven  0 -1  0  0 0 0 0 0;...
        0   0   0 C.cra 0 -C.cra  0  0 -1  0 0 0 0 0;...
        0   0   0   0 C.spi  0    0  0  0 -1 0 0 0 0;...
        0   0   0   0   0    0    1  1  1  0 0 0 0 0;...

        0 0 0 0 0 0 inw 0  0  0 -1  1  0  0;...
        0 0 0 0 0 0  0 inw 0  0  0 -1 -1  0;...
        0 0 0 0 0 0  0  0 inw 0  0  0  0 -1;...
        0 0 0 0 0 0  0  0  0 inw 0  0  0  1;...
        inputA];

    B=[zeros(1,13) inputB]'; % independent terms of the system

    % Constants only affect the mode zero
    if n==0
        B=[ 0 0 -Pjug 0 ...
            -Vo.art -Vo.ven -Vo.cra C.spi*Pout-Vo.spi Vtot ...
            0 0 0 0 inputB]';end

    Xn=A\B; % unknowns of the system, Xn = [P V Q]
    Pn(:,n+1)=Xn(1:6);
    Vn(:,n+1)=Xn(7:10);
    Qn(:,n+1)=Xn(11:end);
end
P=cell2struct(mat2cell(Pn,ones(size(Pn(:,n+1)))'),fieldnames(P),1);
V=cell2struct(mat2cell(Vn,ones(size(Vn(:,n+1)))'),fieldnames(V),1);
Q=cell2struct(mat2cell(Qn,ones(size(Qn(:,n+1)))'),fieldnames(Q),1);
end

function[X]=Reconstruct(Xn,t,w,n)
% This function reconstructs the signals in time domain X from the set of
% Fourier coefficients XN.

fields=fieldnames(Xn);
Xn=struct2cell(Xn);

for i=1:length(Xn) % for each signal
    X.(fields{i})=real(sum(Xn{i}'.*exp(1i*n*w*t)));end % reconstructed signal
end
