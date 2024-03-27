%% BACHELOR'S THESIS - Script 2/3: Model comparison
%
% Thesis title: A compartmental model to investigate intracranial pulsatility
% Author: Josue Perez Sabater (josu.ps13@gmail.com)
% Supervisor: Wilfried Coenen (wcoenen@ing.uc3m.es)
% 
% This script performs the compartmental simulation of intracranial
% pressures, volumes and flows for different model variations (A, B and C).
% Constant parameters and data to reconstruct the signals are defined.
% Matrices are built to represent the system of equations of each model.
% 
% A signal represented as complex Fourier coefficients can be used as an
% additional input to the model. This input signal can be either a flow
% (Qn) measured with MRI or the carotid pressure (Pn).

close all;clear;clc
load ../OUT/Qn.mat % flow coefficients (ml/s)
load ../OUT/Pn.mat % carotid pressure coefficients (mmHg)

%% Define variables and constants
[P,V,Q,Vo,R,C,Pjug,Pout]=Define_VarAndCons(Pn,Qn);

% Data to reconstruct signals
T=3;        % total time (s)
t=0:.01:T;  % time (s)
f=1;        % cardiac frequency (1/s)
w=2*pi*f;   % angular frequency (rad/s)
N=9;        % number of coefficients
n=(0:N-1)'; % coefficient index (0 to 8)
c=colors(); % plot colors

%% Input to the model
input={'Q','art'};

% Solve matrix singularization if needed
singul.activate=1; % set to 1 in order to impose an additional constrain
if singul.activate,singul.input={'Qsas',0};end % average flow through SAS is zero

%% Run models and reconstruct signals
disp('Running model A (full)')
[Pn,Vn,Qn]=RunModel_A(P,V,Q,Vo,R,C,w,N,Pjug,Pout,input);
[p0,v0,q0]=Reconstr(Pn,Vn,Qn,t,w,n,c);

disp('Running model B (simplified 1)')
[Pn,Vn,Qn]=RunModel_B(P,V,Q,Vo,R,C,w,N,Pjug,Pout,input,singul);
[p1,v1,q1]=Reconstr(Pn,Vn,Qn,t,w,n,c);

disp('Running model C (simplified 2)')
[Pn,Vn,Qn]=RunModel_C(P,V,Q,Vo,R,C,w,N,Pjug,Pout,input,singul);
[p2,v2,q2]=Reconstr(Pn,Vn,Qn,t,w,n,c);

set(findobj('type','fig'),'color','w')

%% %%%%%%%%%%%%%%%%%%% FUNCTIONS USED IN THIS SCRIPT %%%%%%%%%%%%%%%%%%%%%%

function[c]=colors()
c.car  =[161   3  24]/255;
c.art  =[208   5   8]/255;
c.arl  =[251 120 184]/255;
c.cap  =[186   3 251]/255;
c.vnl  =[194 120 251]/255;
c.ven  =[  3  90 251]/255;
c.cra  =[  3 119  76]/255;
c.sas  =[ 64 195 143]/255;
c.spi  =[101 253 176]/255;
c.cp_br=[ 96  47  91]/255;
c.cr_br=[ 40  80  44]/255;
c.bra  =[168 168 168]/255;
end

function[P,V,Q,Vo,R,C,Pjug,Pout]=Define_VarAndCons(Pn,Qn)
% Pressures (mm Hg)
P.car=Pn; % pressure in carotid artery
P.art=[]; % pressure in arteries
P.cap=[]; % pressure in capillaries
P.ven=[]; % pressure in veins
P.cra=[]; % pressure in cranial SAS
P.spi=[]; % pressure in spinal SAS
P.bra=[]; % pressure in brain

% Volumes (ml)
V.art=[]; % volume of arteries
V.cap=[]; % volume of capillaries
V.ven=[]; % volume of veins
V.cra=[]; % volume of cranial SAS
V.spi=[]; % volume of spinal SAS
V.bra=[]; % volume of brain

% Flows (ml/s)
Q.art=Qn{3}; % flow through arteries
Q.arl=[];    % flow through arterioles
Q.vnl=[];    % flow through venules
Q.ven=Qn{2}; % flow through veins
Q.sas=Qn{1}; % flow from cranial to spinal SAS
Q.cp_br=[];  % flow from capillaries to brain
Q.cr_br=[];  % flow from cranial SAS to brain

% Compartment volume at rest (ml)
Vo.art=30;   % volume of arteries
Vo.cap=20;   % volume of capillaries
Vo.ven=60;   % volume of veins
Vo.cra=47.5; % volume of cranial SAS
Vo.spi=95;   % volume of spinal SAS
Vo.bra=1400; % volume of brain

% Compartment compliances (ml/mmHg)
C.art=.015; % compliance of arteries
C.cap=.006; % compliance of capillaries
C.ven=.067; % compliance of veins
C.cra=.127; % compliance of cranial SAS
C.spi=.013; % compliance of spinal SAS
C.bra=40;   % compliance of brain

% Vessel resistances (mmHg s/ml)
R.art=1.4;   % resistance in carotid arteries
R.arl=5.6;   % resistance in arterioles
R.vnl=1.7;   % resistance in venules
R.ven=.56;   % resistance in jugular veins
R.sas=.045;  % resistance in cranial-spinal passage
R.cp_br=8e5; % resistance in capillary-brain passage
R.cr_br=800; % resistance in cranial SAS-brain passage

% Other constants
Pjug=3;  % pressure in jugular vein
Pout=10; % pressure out of spinal SAS
end

function[P,V,Q]=RunModel_A(P,V,Q,Vo,R,C,w,N,Pfin,Pout,input)

Vtot=sum(cell2mat(struct2cell(Vo)))-Vo.spi; % total cranial volume

% Find position of input unknown
L=size(input,1);
X=cell(L,1);
PVQ=repmat({P V Q},[L 1]);
[X{:}]=PVQ{cell2mat(input(:,1))=='PVQ'}; % whether the unknown is P, V or Q
unknowns=[strcat('P',fieldnames(P));strcat('V',fieldnames(V));strcat('Q',fieldnames(Q))]';
i=contains(unknowns,string(cell2mat(input))); % position of input

% Run model for each mode
for n=0:N-1 % 0 to 8
    inw=1i*n*w;
    inputA=zeros(L,numel(i));inputA(i)=1; % input coefficients
    inputB=zeros(1,L); % input independent terms
    for x=1:L, inputB(x)=X{x}.(input{x,2})(n+1); end

    % Define system AX=B
    A=[... % coefficients of the system
        1 -1  0  0  0 0  0 0 0 0 0 0 0 -R.art 0   0    0   0    0     0;...
        0  1 -1  0  0 0  0 0 0 0 0 0 0   0 -R.arl 0    0   0    0     0;...
        0  0  1 -1  0 0  0 0 0 0 0 0 0   0    0 -R.vnl 0   0    0     0;...
        0  0  0 -1  0 0  0 0 0 0 0 0 0   0    0   0 -R.ven 0    0     0;...
        0  0  0  0 -1 1  0 0 0 0 0 0 0   0    0   0    0 -R.sas 0     0;...
        0  0  1  0  0 0 -1 0 0 0 0 0 0   0    0   0    0   0 -R.cp_br 0;...
        0  0  0  0  1 0 -1 0 0 0 0 0 0   0    0   0    0   0    0 -R.cr_br;...

        0 C.art 0   0   0   0 -C.art -1  0  0  0  0 0 0 0 0 0 0 0 0;...
        0   0 C.cap 0   0   0 -C.cap  0 -1  0  0  0 0 0 0 0 0 0 0 0;...
        0   0   0 C.ven 0   0 -C.ven  0  0 -1  0  0 0 0 0 0 0 0 0 0;...
        0   0   0   0 C.cra 0 -C.cra  0  0  0 -1  0 0 0 0 0 0 0 0 0;...
        0   0   0   0   0 C.spi  0    0  0  0  0 -1 0 0 0 0 0 0 0 0;...
        0   0   0   0   0   0    0    1  1  1  1  0 1 0 0 0 0 0 0 0;...

        0 0 0 0 0 0 0 inw 0  0  0  0  0 -1  1  0  0  0  0  0;...
        0 0 0 0 0 0 0  0 inw 0  0  0  0  0 -1  1  0  0  1  0;...
        0 0 0 0 0 0 0  0  0 inw 0  0  0  0  0 -1 -1  0  0  0;...
        0 0 0 0 0 0 0  0  0  0 inw 0  0  0  0  0  0 -1  0  1;...
        0 0 0 0 0 0 0  0  0  0  0 inw 0  0  0  0  0  1  0  0;...
        0 0 0 0 0 0 0  0  0  0  0  0 inw 0  0  0  0  0 -1 -1;...
        inputA];

    B=[zeros(1,19) inputB]'; % independent terms of the system

    % Constants only affect the mode zero
    if n==0
        B=[ 0 0 0 -Pfin 0 0 0 ...
            -Vo.art -Vo.cap -Vo.ven -Vo.cra C.spi*Pout-Vo.spi Vtot ...
            0 0 0 0 0 0 inputB]';end

    Xn=A\B; % unknowns of the system, Xn = [P V Q]'
    Pn(:,n+1)=Xn(1:7);
    Vn(:,n+1)=Xn(8:13);
    Qn(:,n+1)=Xn(14:20);
end
P=cell2struct(mat2cell(Pn,ones(size(Pn(:,n+1)))'),fieldnames(P),1);
V=cell2struct(mat2cell(Vn,ones(size(Vn(:,n+1)))'),fieldnames(V),1);
Q=cell2struct(mat2cell(Qn,ones(size(Qn(:,n+1)))'),fieldnames(Q),1);
end

function[P,V,Q]=RunModel_B(P,V,Q,Vo,R,C,w,N,Pfin,Pout,input,singul)

% Remove unused fields
V=rmfield(V,'bra');
Q=rmfield(Q,{'cp_br','cr_br'});
Vo=rmfield(Vo,'bra');
R=rmfield(R,{'cp_br','cr_br'});

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
        1 -1  0  0  0 0  0 0 0 0 0 0 -R.art 0   0    0   0;...
        0  1 -1  0  0 0  0 0 0 0 0 0   0 -R.arl 0    0   0;...
        0  0  1 -1  0 0  0 0 0 0 0 0   0    0 -R.vnl 0   0;...
        0  0  0 -1  0 0  0 0 0 0 0 0   0    0   0 -R.ven 0;...
        0  0  0  0 -1 1  0 0 0 0 0 0   0    0   0    0 -R.sas;...

        0 C.art 0   0   0   0 -C.art -1  0  0  0  0 0 0 0 0 0;...
        0   0 C.cap 0   0   0 -C.cap  0 -1  0  0  0 0 0 0 0 0;...
        0   0   0 C.ven 0   0 -C.ven  0  0 -1  0  0 0 0 0 0 0;...
        0   0   0   0 C.cra 0 -C.cra  0  0  0 -1  0 0 0 0 0 0;...
        0   0   0   0   0 C.spi  0    0  0  0  0 -1 0 0 0 0 0;...
        0   0   0   0   0   0    0    1  1  1  1  0 0 0 0 0 0;...

        0 0 0 0 0 0 0 inw 0  0  0  0 -1  1  0  0  0;...
        0 0 0 0 0 0 0  0 inw 0  0  0  0 -1  1  0  0;...
        0 0 0 0 0 0 0  0  0 inw 0  0  0  0 -1 -1  0;...
        0 0 0 0 0 0 0  0  0  0 inw 0  0  0  0  0 -1;...
        0 0 0 0 0 0 0  0  0  0  0 inw 0  0  0  0  1;...
        inputA];

    B=[zeros(1,16) inputB]'; % independent terms of the system

    % Constants only affect the mode zero
    if n==0
        B=[ 0 0 0 -Pfin 0 ...
            -Vo.art -Vo.cap -Vo.ven -Vo.cra C.spi*Pout-Vo.spi Vtot ...
            0 0 0 0 0 inputB]';end

    Xn=A\B; % unknowns of the system, Xn = [P V Q]'
    Pn(:,n+1)=Xn(1:7);
    Vn(:,n+1)=Xn(8:12);
    Qn(:,n+1)=Xn(13:end);
end
P=cell2struct(mat2cell(Pn,ones(size(Pn(:,n+1)))'),fieldnames(P),1);
V=cell2struct(mat2cell(Vn,ones(size(Vn(:,n+1)))'),fieldnames(V),1);
Q=cell2struct(mat2cell(Qn,ones(size(Qn(:,n+1)))'),fieldnames(Q),1);
end

function[P,V,Q]=RunModel_C(P,V,Q,Vo,R,C,w,N,Pjug,Pout,input,singul)

% Remove unused fields
P=rmfield(P,'cap');
V=rmfield(V,{'cap','bra'});
Q.cap=[];
Q=rmfield(Q,{'arl','vnl','cp_br','cr_br'});
Q=orderfields(Q,[1 4 2 3]);
Vo=rmfield(Vo,{'cap','bra'});
R.cap=R.arl+R.vnl;
R=rmfield(R,{'arl','vnl','cp_br','cr_br'});

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

function[P,V,Q]=Reconstr(Pn,Vn,Qn,t,w,n,c)
% This function reconstructs the pressure, volume and flow signals in time
% domain from the set of Fourier coefficients of each one. All signals are
% plotted against time.

Xn={Pn Vn Qn}; % signals as Fourier coefficients
X=cell(size(Xn)); % signals in time domain
Xn_names={'Pressure (mmHg)','Volume (cm^3)','Flow (cm^3/s)'};
Xn_titles={'Pressure waveform','Volume waveform','Flow waveform'};

figure
for x=1:length(Xn) % for each magnitude
    subplot(1,3,x),hold on
    fields=fieldnames(Xn{x});
    Xn{x}=struct2cell(Xn{x});

    for i=1:length(Xn{x}) % for each signal
        X{x}.(fields{i})=real(sum(Xn{x}{i}'.*exp(1i*n*w*t))); % reconstructed signal
        plot(t,X{x}.(fields{i}),'Color',c.(fields{i}))
        xlabel('Time (s)'),end
    title(Xn_titles{x}),ylabel(Xn_names{x}),legend(fields)

    if x==2,set(gca,'YScale','log'),end % volume is plotted in log scale
end
P=X{1};V=X{2};Q=X{3};
end
