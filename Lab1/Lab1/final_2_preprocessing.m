%% S. Rakotonarivo
% Created on July 24th, 2025
% Last modified Nov 13th, 2025
% Validated on Nov 13th, 2025
%
% This code:
%  - loads data in .mat format , 
%  - applies a bandpass filter on them,
%  - applies the welch algorithm + Xcorr for predicting Green's functions
%  - stuctures the result in a 3D matrix that contains the Green's function for
% all pairs of mics in time domain,
%  - and saves the resulting matrix as well as other parameter in a structure "I" in  a mat file.
%
% The structure as the following fields:
% I.C_t         2Nbmic x 2Nbmic x ns     Xcorrelation
% I.Gt          2Nbmic x 2Nbmic x ns     Green's function 
% I.t_Gt        1xns    time vector for Green's function
% I.fs          1x1     sampling frequency
% I.Nbmic       1x1     Number of pairs of mics 
% I.ns          1x1     Number of samples
% I.duration    1x1     Duration Green's function
% I.n_segment   1x1     Number of time segment used for computing Xcorrelation
% I.c0          1x1     Sound speed during acquisition

clearvars; %close all
%------------------------------
%% Initialize variables.
%------------------------------

% --- PATH SETTINGS (Student Updated) ---
Generaldir = 'D:\france semester\3_Acoustic Imaging\lab\Shahariar\DATA\';
FileDir = strcat(Generaldir, 'Data_MAT\');

% --- CONFIGURATION (Select One) ---
% FileDataName = 'MES_RealCylinder_';   FileSaveName = 'CC50_Cylinder.mat'; % Config 1
FileDataName = 'MES_RealEmpty_';      FileSaveName = 'CC50_Empty.mat';    % Config 2

FlagSave = 1; 

% Parameters: Acquisition 
Nbmic = 16;       % Number of pairs of microphones (Total 32)
nReal = 50;       % Number of files to load
ii0 = 1;          % Index of first test to load

% Sound speed
Temperature = 23;
c0 = 343;

% Filter parameters
fmin = 1000; % Low Frequency [Hz]
fmax = 9500; % High Frequency [Hz]

% Parameter for Welch algorithm 
ns_short = 1638;    % Length of sub-signal
poverlap = 80;      % Overlap percent

%% ----------------------------------------
% Get parameters and Design Filter
%-----------------------------------------

% Load first data file to get sampling frequency
ii = 1;
FileName = strcat(FileDir, FileDataName, num2str(ii), '.mat');
if ~isfile(FileName), error('File not found. Check LoadData step.'); end
load(FileName); % Loads variable 'A'

% Parameters
fs = A.fsample; 
ns = length(A.t);
df = fs/ns;
dt = 1/fs;
t = A.t; 

% Bandpass filter design (FIR)
disp('Designing FIR Bandpass Filter...');
filter = designfilt('bandpassfir', 'FilterOrder', 2000, ... 
     'CutoffFrequency1', round(fmin/df)*df, 'CutoffFrequency2', round(fmax/df)*df, ...  
     'SampleRate', fs);

% Data matrix initialization
data_filt = zeros(Nbmic*2, nReal, ns); % 32 mics x 50 realizations x Samples

%%----------------------------------------
% Data loading & filtering & structuration: to be coded by STUDENTS
%-----------------------------------------
disp('Filtering Data...');

for ii = ii0:1:ii0+nReal-1  % Loop over each data acquisition
    tic
        % ---- 1) Load data ----
        FileName = strcat(FileDir, FileDataName, num2str(ii), '.mat');
        load(FileName);
        
        % Concatenate Left (Outer) and Right (Inner) microphones
        % A.yl is 16xN, A.yr is 16xN -> Result is 32xN
        raw_signal = [A.yl; A.yr];

        % ---- 2) Filter data ----
        % Loop through all 32 channels
        for mic = 1:(Nbmic*2)
            % Apply zero-phase filtering to align signals perfectly
            data_filt(mic, ii, :) = filtfilt(filter, raw_signal(mic, :));
        end
   
     % sprintf('Acquisition n°%s/%s: filtering performed in %f s',num2str(ii),num2str(nReal),toc)
end
disp('Filtering Complete.');

%% ----------------------------------------
% Average Xcorrelation --> Green's functions estimation
%-----------------------------------------

% Average Xcorrelation calculation
disp('Calculating Cross-Correlation (Welch Method)...');
tic
% Using the student version to ensure Real outputs (fixes complex noise)
[C_t, duration_segment, n_segment] = FCT_Xcorr_Welch_Student(data_filt, ns_short, fs, poverlap);
toc

% Matrix conditioning for estimating the Green's functions
% Applying Equation (2) from the manual (Time Derivative)
dt_Xcorr = 1/fs;
dC_t = diff(C_t, [], 3);
Gt = dC_t(:, :, ns_short/2+1:end) ./ dt; % Keep positive time (causal part)
t_Gt = 0:1/fs:(length(Gt)-1)/fs;

%% ----------------------------------------
% Data saving
%-----------------------------------------
if FlagSave == 1
    I.C_t = C_t;     % Xcorrelation
    I.Gt = Gt;       % Green's function 
    I.t_Gt = t_Gt;   % Time vector
    I.fs = fs;
    I.Nbmic = Nbmic;
    I.ns = length(t_Gt);
    I.duration = (duration_segment-dt*2)/2;
    I.n_segment = n_segment;
    I.c0 = c0;
    
    SaveName = strcat(FileDir, FileSaveName);
    save(SaveName, 'I', '-v7.3');
    disp(['Saved processed data to: ', SaveName]);
end


%% HELP ON function FCT_Xcorr_Welch:
% 
% [Ct,duration_short,nb_segment] = FCT_Xcorr_Welch(data_t,ns_short,fs,poverlap)
%
% S Rakotonarivo, July 31st 2025
% Last modified July 31st 2025
% Validated July 31st 2025   
%
%******************************************************
% LABORATORY WORK: WELCH METHOD (portion of codes)
%******************************************************
%
%------------------------
% INPUT PARAMETERS
%------------------------
% data_t     nMic x nReal x nSamples     time domain signal 
% poverlap   1 x 1      natural percent overlapping between window for short Fourier transforms
% ns_short   1 x 1      number os sample of the signal sub-portion/sub-segment
% fs         1 x 1      sampling frequency (Hz)
%------------------------
% OUTPUT PARAMETERS
%------------------------
% Ct         nMic x nMic x nSamplesSegment    Normalized Noise correlation function
% duration_short   1 x 1      time duration of sub-segment /sub-portion of signal
% nb_segment       1 x 1      number of sub-segments/sub-portions of signal 



