%% BACHELOR'S THESIS - Script 1/3: Data processing
%
% Thesis title: A compartmental model to investigate intracranial pulsatility
% Author: Josue Perez Sabater (josu.ps13@gmail.com)
% Supervisor: Wilfried Coenen (wcoenen@ing.uc3m.es)
%
% The data for this script, stored in the folder "DCM", consists of 240 MRI
% images taken from the same patient at different velocity encodings (VENC)
% organized as follows:
% - 4 folders corresponding to different series (mainly different VENC)
% - 60 images in each folder (30 magnitude and 30 phase) that together
%   correspond to 30 time instants of a reconstructed cardiac cycle
%
% Three series contain information of the same body plane (C2 vertebrae),
% while the fourth series contains information of the aqueduct. Four
% regions are defined, which provide information on blood and CSF flow: CSF
% (CSF through foramen magnum), VEN (venous blood), ART (arterial blood)
% and CER (CSF through cerebral aqueduct).
% - C2-V10: CSF, VEN and ART regions
% - C2-V30: CSF, VEN and ART regions
% - C2-V80: CSF, VEN and ART regions
% - AQ-V15: CER region
%
% Each region may contain several vessels to study. Variables with data of
% different vessels will be organized in matrices M{r,v} where "r" is the
% region and "v" the vessel. Additionally, each matrix element may contain
% data for different series "s" and time instants "t": M{r,v}{s}{t}
%
% The script returns a variable "Q" with the flow (cm^3/s) through each
% region along one cardiac cycle, and a variable "t" containing the time.

close all;clear;clc

% Define directories and data organization
dcmdir='../DCM/'; % contains DICOM files
outdir='../OUT/'; % contains output (results)
viddir='../VID/'; % contains output (videos)
SER={'C2-V10','C2-V30','C2-V80','AQ-V15'}; % series (LO ME HI AQ)
REG={'CSF';'VEN';'ART';'CER'};             % regions
VES={'SSAS','DURA','CORD','    ';          % vessels
     'LIJV','RIJV','LEDV','REDV';
     'LICA','RICA','LVEA','RVEA';
     'AQSY','    ','    ','    '};
pwdlist={[dcmdir SER{1}]... % pathway for each series
    [dcmdir SER{2}]...
    [dcmdir SER{3}]...
    [dcmdir SER{4}]};

nt10=define_colors; % colors used in some plots

% Set to 1 in order to perform the action:
opt.reaDCM=0; % read data from DICOM files (otherwise load)
opt.wriMOV=0; % write movies with velocity field
opt.defROI=0; % manually define regions of interest (otherwise load)
opt.shoROI=0; % show regions of interest
opt.repSIG=0; % repeat the signal
opt.pltSER=0; % plot flow of all series
opt.savALL=0; % save all workspace
opt.savFLW=0; % save flow signals

if opt.repSIG,rep=4; % number of times to repeat signal
else,rep=1;end

% Define Fourier transform parameters
f=1;           % cardiac frequency (1/s)
w=2*pi*f;      % angular frequency (rad/s)
N=9;           % number of Fourier coefficients
t_rec=0:.01:1; % time for reconstructed signal (s)

% Structure I contains indeces of different levels:
% I.ser contains the series index:                      (1 2 3 4)
%
% I.reg determines which regions are analyzed in each series:
% - Series C2-V10 contains regions CSF, VEN, ART        (1 2 3  )
% - Series C2-V30 contains regions CSF, VEN, ART        (1 2 3  )
% - Series C2-V80 contains regions CSF, VEN, ART        (1 2 3  )
% - Series AQ-V15 contains region CER                   (      4)
%
% I.ves determines which vessels are analyzed in each region:
% - Region CSF contains vessels SSAS, DURA, CORD        (1 2 3  )
% - Region VEN contains vessels LIJV, RIJV, LEDV, REDV  (1 2 3 4)
% - Region ART contains vessels LICA, RICA, LVEA, RVEA  (1 2 3 4)
% - Region CER contains vessel  AQSY                    (1      )
I.ser=1:4;             % index of series
I.reg={1:3;1:3;1:3;4}; % index of regions analyzed in each series
I.ves={1:3;1:4;1:4;1}; % index of vessels analyzed in each region

%% Read (or load) DICOMs, write movies and reconstruct magnitude images
if opt.reaDCM
    DAT=read_dicoms(I,pwdlist,outdir);
else;load ../OUT/DAT.mat;end
I.tim=1:DAT.Nt; % time index

if opt.wriMOV
    make_movies(SER,I,DAT.venc,DAT.time,DAT.U,viddir);end
if opt.defROI||opt.shoROI
    IMG=reconstruct_IMG(SER,I,DAT.magni);end

%% Define (or load) ROIs in each series, show them and combine them
if opt.defROI
    ROI=define_ROI(I,IMG,outdir);
else;load ../OUT/ROI.mat;end

if opt.shoROI
    show_ROI(SER,REG,I,IMG,ROI);end
ROI_C=combine_ROI(I,ROI);

%% Get properties of ROIs
[PRO,PRO_C]=get_props_ROI(I,ROI,ROI_C);
[A,U,U_avg]=apply_ROI(I,DAT.PXA,DAT.U,ROI);

%% Correct velocity and compute (and repeat) flow Q of each ROI
U_corrt=remove_alias(I,DAT.venc,ROI_C,U,PRO,PRO_C);
[time,Q,shift]=compute_Q(I,DAT.time,DAT.T,DAT.PXA,U_corrt);
if opt.repSIG==1
    [time,Q]=repeat_Q(I,time,Q,rep);end

%% Plot (and save) the results
if opt.pltSER==1
    plot_all_SER(time,DAT.PXA,U,Q,rep,shift,SER,REG,VES);end
Q=plot_Q(time,Q,nt10);

%% Represent signal as Fourier coefficients
t=time(I.tim)/100/f;           % time of original signal
t=t-t(1);                      % shift time to start at 0
for r=1:3;Q{r}=Q{r}(I.tim);end % remove repeated signal, take only one CC
Qn=FourierCoeff(REG,t,Q,N,w,t_rec);

% Carotid pressure waveform
a=[1 -0.0345 -0.0511 -0.0267 -0.0111 -0.0013 0.0050 0.0027  0.0061]';
b=[0 -0.1009 -0.0284  0.0160  0.0070  0.0174 0.0041 0.0041 -0.0005]';
Pn=a+b*1i; Pn=102.3530*Pn;

%% Save and print results
set(findobj('type','fig'),'color','w')
if opt.savALL
    save(outdir+"results.mat");end
if opt.savFLW
    save(outdir+"Qn.mat",'Qn');save(outdir+"Pn.mat",'Pn');end

fprintf("\n")
disp("Average arterial flow: "+mean(Q{3})+" cm^3/s")
disp("Average venous flow: "  +mean(Q{2})+" cm^3/s")
disp("Average CSF flow: "     +mean(Q{1})+" cm^3/s")
disp("Net arterial flow: " +trapz(Q{3})+" cm^3/CC")
disp("Net venous flow: "   +trapz(Q{2})+" cm^3/CC")
disp("Net CSF flow: "      +trapz(Q{1})+" cm^3/CC")

%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS USED IN THIS SCRIPT %%%%%%%%%%%%%%%%%%%%%%

%% Define colors
function nt10=define_colors()
% Colors to be used in the plots

nt10.blue   = [0.30588 0.47451 0.65490];
nt10.green  = [0.34902 0.63137 0.30980];
nt10.brown  = [0.61176 0.45882 0.37255];
nt10.orange = [0.94902 0.55686 0.16863];
nt10.yellow = [0.92941 0.78824 0.28235];
nt10.gray   = [0.72941 0.69020 0.67451];
nt10.red    = [0.88235 0.34118 0.34902];
nt10.purple = [0.69020 0.47843 0.63137];
nt10.teal   = [0.46275 0.71765 0.69804];
nt10.pink   = [1.00000 0.61569 0.65490];
end

%% Read DICOMs
function[DAT]=read_dicoms(I,pwdlist,outdir)
% Get and store data in cell arrays where "s" is the series and "i"/"t" is
% the image/time instant. Each series contains 60 images (30 magnitude and
% 30 phase) which are sorted first all magnitude and then all phase.
% Variables marked with an asterisk (*) are stored in the structure DAT.
disp("Reading DICOM images...")

S=length(I.ser);         % number of series

info         =cell(1,S); % all data of each file
img          =cell(1,S); % all images (60)
phase        =cell(1,S); % phase images
magni        =cell(1,S); % magnitude images                      *
U            =cell(1,S); % velocity field images (cm/s)          *
flag         =cell(1,S); % magnitude/phase flag (8/11)
trigger_time =cell(1,S); % acquisition time (ms)
time         =cell(1,S); % acquisition time (% of CC)            *

T           =zeros(1,S); % length of cardiac cycle (ms)          *
Nt          =zeros(1,S); % number of acquisition time instants   *
scale       =zeros(1,S); % velocity encode scale
venc        =zeros(1,S); % velocity encoding (cm/s)              *
venc_S      =zeros(1,S); % velocity encoding (scaled to decode phase image)
loc         =zeros(1,S); % slice location (mm)                   *
PXA         =zeros(1,S); % pixel area (cm^2)                     *

for s=I.ser
    dicomlist=dir(fullfile(pwdlist{s},'*.dcm')); % all DICOM images in series S
    N=length(dicomlist); % number of images in series S

    % Collect data from files and sort it by flag (first magnitude, then
    % phase) and then by triggertime (ascending)
    for i=1:N
        info{s}{i}=dicominfo(fullfile(pwdlist{s},dicomlist(i).name));
        img{s}{i}=double(dicomread(info{s}{i}));
        flag{s}(i)=info{s}{i}.(dicomlookup('0043','1030'));
        trigger_time{s}(i)=info{s}{i}.TriggerTime;end % 0018,1060

    [~,j]=sortrows([flag{s}',trigger_time{s}'],[1 2]);
    img{s}=img{s}(j);

    % Collect data unique for each series, but common for all time instants
    T(s)      =info{s}{1}.NominalInterval; % 0018,1062
    time{s}   =100*unique(trigger_time{s}(j))/T(s);
    Nt(s)     =length(time{s});
    scale(s)  =double(info{s}{1}.(dicomlookup('0019','10E2')));
    venc(s)   =double(info{s}{1}.(dicomlookup('0019','10CC')));
    venc_S(s) =pi*scale(s)/venc(s);
    loc(s)    =info{s}{1}.SliceLocation; % 0020,1041
    PXA(s)    =prod(info{s}{1}.PixelSpacing)/100; % 0028,0030

    % Convert the phase and magnitude images into velocity field
    for t=1:Nt(s)
        magni{s}{t}=img{s}{t};
        phase{s}{t}=img{s}{t+Nt(s)};
        U{s}{t}=-.1*(phase{s}{t}./max(magni{s}{t},1))/venc_S(s);end;end

% Save relevant variables. Time data (time, T, Nt) is assumed to be common
% for all series, and so it is taken from first series (C2-V10)
DAT.time=time{1};    DAT.venc=.1*venc;       DAT.PXA=PXA;
DAT.T=T(1);          DAT.magni=magni;        DAT.loc=loc;
DAT.Nt=Nt(1);        DAT.U=U;

save([outdir,'DAT.mat'],'DAT')
end

%% Make movies
function[]=make_movies(SER,I,venc,time,U,viddir)
disp("Making movies...")
fg=figure;
for s=I.ser
    vid=VideoWriter([viddir 'U_tot_' SER{s}]);
    vid.Quality=100;
    vid.FrameRate=20;
    open(vid);

    for t=I.tim
        figure(fg)
        pcolor(U{s}{t})
        shading flat
        axis equal
        set(gca,'YDir','reverse')
        clim(venc(s)*[-1 1]);colorbar
        title("Series "+SER{s})
        xlabel("t = "+num2str(time(t),'%.2f')+"% of CC")
        drawnow
        pause(.5)

        frame=getframe(fg);
        writeVideo(vid,frame);end
    close(vid);end
end

%% Reconstruct magnitude images
function IMG=reconstruct_IMG(SER,I,magni)
% IMG contains a magnitude image for each series, reconstructed from 30
% magnitude MRI files. First, we define IMG by taking the highest value at
% each pixel position among all 30 images. Then, we normalize IMG (dividing
% it by its maximum value) and map the values to new intensities: 0 to 0.8
% maps to 0 to 1 (high-intensity values become saturated, increasing
% contrast in low-intensity values).
%
% The figures displayed can be used to compare the position of vessels in
% each series, and check if there is displacement. It is also used to
% manually select the ROIs, in case that option is activated.

disp("Reconstructing magnitude images...")
IMG=cell(size(I.ser));

% Reconstruct and display magnitude image.
fg(1)=figure;
for s=I.ser
    IMG{s}=0;
    for t=I.tim
        IMG{s}=max(IMG{s},magni{s}{t});end
    IMG{s}=imadjust(IMG{s}/max(IMG{s}(:)),[0 .8]);

    subplot(2,2,s)
    imshow(IMG{s})
    title("Series "+s+": "+SER{s});end

% Display a combination of magnitude images from different series, useful
% to detect displacement of vessels between series.
fg(2)=figure;
subplot(2,2,1)
imshow(imfuse(IMG{1},IMG{2}))
title("Series 1 and 2: "+SER{1}+" and "+SER{2})
subplot(2,2,2)
imshow(imfuse(IMG{2},IMG{3}))
title("Series 2 and 3: "+SER{2}+" and "+SER{3})
subplot(2,2,3)
imshow(imfuse(IMG{3},IMG{1}))
title("Series 3 and 1: "+SER{3}+" and "+SER{1})
subplot(2,2,4)
imshow(cat(3,IMG{1},IMG{2},IMG{3}))
title("Series 1, 2 and 3")

linkaxes(findall(fg,'type','axes'))
end

%% Define regions of interest (ROI)
function ROI=define_ROI(I,IMG,outdir)
% ROIs are manually selected on the screen, on top of the reconstructed
% magnitude image IMG. Here, the "for" loop that goes through each series
% must be the last one, so indeces are not the same as in variable I:
I.reg=1:4;             % index of regions
I.ves={2:3;1:4;1:4;1}; % index of vessels analyzed in each region
I.ser={1:3;1:3;1:3;4}; % index of series containing each region

disp("Defining ROIs")
VES_name={'SSAS'...
    'DURA MATER'...
    'PIA MATER','';
    'LEFT INTERNAL JUGULAR VEIN'...
    'RIGHT INTERNAL JUGULAR VEIN'...
    'LEFT EPIDURAL VEIN'...
    'RIGHT EPIDURAL VEIN';
    'LEFT INTERNAL CAROTID ARTERY'...
    'RIGHT INTERNAL CAROTID ARTERY'...
    'LEFT VERTEBRAL ARTERY'...
    'RIGHT VERTEBRAL ARTERY';
    'CEREBRAL AQUEDUCTS','','',''};
str={'all 3 images';'all 3 images';'all 3 images';'the image'};
ROI=cell(size(VES_name));

% Draw the ROI of each vessel (except SSAS)
figure,fg=get(gcf,'Number');
for r=I.reg
    for v=I.ves{r}
        disp("Please contour the "+VES_name{r,v}+" on "+str{r})
        for s=I.ser{r}
            figure(fg)
            imshow(IMG{s},'InitialMagnification',400);colorbar
            ROI{r,v}{s}=createMask(drawpolygon());
            close(fg);end;end;end

% Compute ROI of vessel SSAS (DURA - CORD)
for s=I.ser{1};ROI{1,1}{s}=ROI{1,2}{s}-ROI{1,3}{s};end

save([outdir,'ROI.mat'],'ROI')
end

%% Show ROI images
function[]=show_ROI(SER,REG,I,IMG,ROI)
% Display ROIs over the combined magnitude images IMG

disp('Displaying ROIs...')
I.ves{1}=1; % in CSF region, only show SSAS vessel

fg=figure().Number;
for s=I.ser
    for r=I.reg{s}
        ROI_REG=0; % ROI of all vessels in the same region combined
        for v=I.ves{r};ROI_REG=ROI_REG+ROI{r,v}{s};end

        figure(fg+s-1);subplot(2,2,r)
        imshow(imfuse(IMG{s},ROI_REG))
        title(REG{r});end
    set(gcf,'Name',"ROIs in series "+SER{s});end
end

%% Combine ROI
function[ROI_C]=combine_ROI(I,ROI)
% ROI_C{R,V} is the combination of ROIs of all series in region R, vessel
% V. There are three methods to combine ROIs: obtaining the union of all
% series, their intersection, or just choosing a specific series. Comment
% all but the option to use

% % Combine ROI: Union
% ROI_C=num2cell(zeros(size(ROI)));
% for s=I.ser
%     for r=I.reg{s}
%         for v=I.ves{r}
%             ROI_C{r,v}=(ROI_C{r,v}+ROI{r,v}{s})>0;end;end;end

% Combine ROI: Intersection
ROI_C=num2cell(ones(size(ROI)));
for s=I.ser
    for r=I.reg{s}
        for v=I.ves{r}
            ROI_C{r,v}=ROI_C{r,v}.*ROI{r,v}{s};end;end;end

% % Combine ROI: Use ROI of one series
% s=[1 1 1 4]; % choose series to represent the combined ROI for each region
% for r=1:4
%     for v=I.ves{r}
%         ROI_C{r,v}=ROI{r,v}{s(r)};end;end
end

%% Get properties from ROI
function[PRO,PRO_C]=get_props_ROI(I,ROI,ROI_C)
% Obtain properties stated in PROPS from all the ROIs, including the
% combined ROI (ROI_C). Properties are returned in PRO and PRO_C.

props={'Area','BoundingBox','Centroid','EquivDiameter','MajorAxisLength','MinorAxisLength','Orientation'};
PRO=cell(size(ROI));
PRO_C=cell(size(ROI_C));

for s=I.ser
    for r=I.reg{s}
        for v=I.ves{r}
            PRO{r,v}{s}=regionprops(ROI{r,v}{s},props);
            cent=round(PRO{r,v}{s}.Centroid);
            PRO{r,v}{s}.x=cent(1);
            PRO{r,v}{s}.y=cent(2);end;end;end

for r=1:4
    for v=I.ves{r}
        PRO_C{r,v}=regionprops(ROI_C{r,v},props);
        cent=round(PRO_C{r,v}.Centroid);
        PRO_C{r,v}.x=cent(1);
        PRO_C{r,v}.y=cent(2);end;end
end

%% Get area and velocity of ROI
function[A,U_roi,U_avg]=apply_ROI(I,PXA,U,ROI)
% Obtain more properties of ROIs: area, velocity and mean velocity.

A=cell(size(ROI));     % area (cm^2) in each vessel and series
U_roi=cell(size(ROI)); % velocity in each ROI (vessels), series and time instant
U_avg=cell(size(ROI)); % average U in each ROI, series and time instant

for s=I.ser
    for r=I.reg{s}
        for v=I.ves{r}
            A{r,v}{s}=sum(ROI{r,v}{s}(:))*PXA(s);
            for t=I.tim
                U_roi{r,v}{s}{t}=U{s}{t}.*ROI{r,v}{s};
                U_avg{r,v}{s}(t)=mean(nonzeros(U_roi{r,v}{s}{t}));end;end;end;end
end

%% Remove aliasing in U
function[U_corrt]=remove_alias(I,venc,ROI_C,U,PRO,PRO_C)
% U_CORRT is a new velocity field where aliasing has been removed. Data
% from all series is combined to determine the true velocity at each pixel
% within the combined ROI (ROI_C). U_CORRT is not computed for region CER.

disp("Removing alising...")
I.reg=1:3;             % index of regions
I.ves={1:3;1:4;1:4};   % index of vessels analyzed in each region
I.ser={1;1:3;1:3};     % index of series analyzed in each region
U_trans=cell(size(U)); % translated velocity field
U_corrt=cell(size(U)); % corrected velocity field
vLO=-8*venc(1):2*venc(1):8*venc(1);
vME=-4*venc(2):2*venc(2):4*venc(2);

% Remove aliasing - CSF
for r=I.reg(1)
    for v=I.ves{r}
        for s=I.ser{r}
            for t=I.tim
                U_corrt{r,v}{t}=U{r,v}{s}{t};end;end;end;end

% Remove aliasing - VEN and ART
for r=I.reg(2:3)
    for v=I.ves{r}
        [row,col]=find(ROI_C{r,v}); % coordinates of ROI_C pixels
        for s=I.ser{r} % translation vector from each series ROI to ROI_C
            trans{s}=[PRO_C{r,v}.x-PRO{r,v}{s}.x,PRO_C{r,v}.y-PRO{r,v}{s}.y];end

        for t=I.tim
            U_corrt{r,v}{t}=zeros(size(U{r,v}{1}{t}));
            for s=I.ser{r} % velocity field of translated ROI
                U_trans{r,v}{s}{t}=imtranslate(U{r,v}{s}{t},trans{s});end

            for k=1:length(row)
                rk=row(k);
                ck=col(k);

                [~,ind]=min(abs(U_trans{r,v}{2}{t}(rk,ck)-U_trans{r,v}{3}{t}(rk,ck)+vME));
                U_trans{r,v}{2}{t}(rk,ck)=U_trans{r,v}{2}{t}(rk,ck)+vME(ind);

                [~,ind]=min(abs(U_trans{r,v}{1}{t}(rk,ck)-U_trans{r,v}{2}{t}(rk,ck)+vLO));
                U_corrt{r,v}{t}(rk,ck)=U_trans{r,v}{1}{t}(rk,ck)+vLO(ind);end;end;end;end

% For region VEN, vessel LEDV, we take velocity in series C2-V10 as the
% corrected velocity
U_corrt{2,3}=U_trans{2,3}{1};
end

%% Compute flow Q of each ROI
function[time,Q,j]=compute_Q(I,time,T,PXA,U)
% Compute flow Q in each vessel and shift time reference. Q is not computed
% for series AQ-V15 (region CER). Time in DICOM files is referenced to the
% finger ECG. We shift the time variable by summing the heart-to-finger PTT
% to reference it to the heart ECG.

% PTT is scaled and added to TIME, and points exceeding 100% are moved to
% the beginning of the vector
deltaT=220; % pulse transit time (PTT) from heart to finger (ms)
time=time+100*deltaT/T;
time(time>100)=time(time>100)-100;
[time,j]=sort(time);

Q=cell(size(U)); % flow (cm^3/s)
for r=1:3
    for v=I.ves{r}
        for t=I.tim
            Q{r,v}(t)=sum(U{r,v}{t}(:))*PXA(1);end
        Q{r,v}=Q{r,v}(j);end;end
end

%% Repeat the flow
function[time_rep,Q_rep]=repeat_Q(I,time,Q,rep)
% Repeat the Q signal REP times for better visualization, as well as the
% time dimension.

disp("Repeating the signal a few times...")
Nt=numel(time); % number of time instants
Q_rep=cell(size(Q)); % flow repeated REP times
time_rep=repmat(time,1,rep)+repelem(0:rep-1,Nt)*100; % time extended REP times

for r=1:3
    for v=I.ves{r}
        Q_rep{r,v}=repmat(Q{r,v},1,rep);end;end
end

%% Plot flow of all series and corrected flow
function[]=plot_all_SER(time,PXA,U,Q,rep,shift,SER,REG,VES)
% Compute, shift, repeat and plot flow through each vessel with each
% velocity encoding (Q_VENC), together with the corrected flow (Q) with no
% aliasing.

disp("Plotting flows with each velocity encoding...")
I.reg=1:3;           % index of regions
I.ves={1;1:4;1:4};   % index of vessels analyzed in each region
I.ser={1:3;1:3;1:3}; % index of series analyzed in each region
I.tim=1:30;          % time index

Q_venc=cell(size(U)); % flow with each VENC (cm^3/s)
cmp=[                 % blue colormap
    .8706 .9216 .9686;  % light
    .6196 .7922 .8824;  % middle
    .1922 .5098 .7412]; % dark

fg=figure().Number;
for r=I.reg
    figure(fg+r-1)
    for v=I.ves{r}
        subplot(length(I.ves{r}),1,v),hold on
        title(REG{r}+" "+VES{r,v})

        for s=I.ser{r}
            for t=I.tim
                Q_venc{r,v}{s}(t)=sum(U{r,v}{s}{t}(:))*PXA(s);end % compute Q
            Q_venc{r,v}{s}=Q_venc{r,v}{s}(shift);                 % shift Q
            Q_venc{r,v}{s}=repmat(Q_venc{r,v}{s},1,rep);          % repeat Q
            plot(time,Q_venc{r,v}{s},'.-','Color',cmp(s,:));end   % plot Q

        plot(time,Q{r,v},'r-')
        legend({SER{I.ser{r}},'Corrected Q'});end;end
end

%% Plot the flow
function[Q_T]=plot_Q(time,Q,nt10)
% Q_T contains the total flow of each region, at each time instant. Here
% the first "for" loop goes through all regions, so indeces are not the
% same as in variable I:

disp("Plotting total flows...")
I.reg=1:3;         % index of regions
I.ves={1;1:4;1:4}; % index of vessels analyzed in each region

% Compute total flow of each region
Q_T=cell(size(I.reg)); % total flow (cm^3/s)
for r=I.reg
    Q_T{r}=0;
    for v=I.ves{r}
        Q_T{r}=Q_T{r}+Q{r,v};end;end
Q_rem=-(Q_T{1}+Q_T{2}+Q_T{3}); % remainder venous flow
Q_T{2}=Q_T{2}+Q_rem;

% Plot total flow of each region
figure
title('Flow into the cranial vault (II)')
hold on;grid on
yline(0,'--')
plot(time,Q_T{1},'Color',nt10.teal,'linewidth',2)   % CSF flow
plot(time,Q_T{3},'Color',nt10.red,'linewidth',2)    % arterial flow
plot(time,Q_T{2},'Color',nt10.purple,'linewidth',2) % venous flow
plot(time,Q_T{2}-Q_rem,'Color',[nt10.purple .5]);   % jugular venous flow
plot(time,Q_rem,'--','Color',[nt10.purple .5]);     % remainder venous flow
legend('','CSF','A','V','V_j','V_r')
end

%% Transform signal to Fourier coefficients
function[Xn]=FourierCoeff(REG,t,X,N,w,t_rec)
% FOURIERCOEFF returns a column XN{r} with the N first coefficients of the
% Fourier series that is equivalent to the signal X{r}. The number of
% Fourier coefficients must be equal or smaller than the length of the
% original signal.

disp("Decomposing signal in "+N+" Fourier coefficients...")
Xn=cell(size(X));  % signal coefficients

figure
for r=1:3
    subplot(1,3,r),hold on
    N=min(N,length(X{r}));     % number of coefficients
    n=(0:N-1)';                % coefficient index (C0, C1, ...)
    Y=fft(X{r})'/length(X{r}); % all Fourier coefficients
    Xn{r}=Y(n+1);              % first N coefficients
    Xn{r}(2:end)=2*conj(Xn{r}(2:end));

    X_rec=real(sum(Xn{r}.*exp(1i*n*w*t_rec))); % reconstructed signal

    plot(t,X{r},t_rec,X_rec)
    title("Region "+REG{r})
    xlabel('Time (% of CC)'),ylabel('Flow (cm^3/s)')
    legend('Original signal','Reconstructed signal');end
end
