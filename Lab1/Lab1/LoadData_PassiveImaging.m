%% S. Rakotonarivo
% Modified for Student Lab - Acoustic Imaging
% This code loads I2S data and saves them in .mat format.

clearvars; close all

%------------------------------------------------
%% Parameters to be changed/updated by students
%-------------------------------------------------

% --- PATH SETTINGS (Updated for your Computer) ---

% Base Directory containing your Config folders
Generaldir = 'D:\france semester\3_Acoustic Imaging\lab\Shahariar\DATA\';

% Output Directory (Where .mat files will be saved)
% IMPORTANT: Make sure you create a folder named "Data_MAT" inside DATA first!
SaveDir = strcat(Generaldir, 'Data_MAT\');

% -----------------------------------------------------------
% CHOOSE CONFIGURATION (Comment/Uncomment lines below)
% -----------------------------------------------------------

% OPTION A: RUN THIS FOR CONFIG 1 (WITH CYLINDER)
% FileDir = strcat(Generaldir, 'Config1\');     % Folder for Config 1
% FileSaveName = 'MES_RealCylinder_';           % Output filename for Cylinder

% OPTION B: RUN THIS FOR CONFIG 2 (WITHOUT CYLINDER)
% (After finishing Config 1, uncomment these two lines and comment out Option A)
FileDir = strcat(Generaldir, 'Config2\');   % Folder for Config 2
FileSaveName = 'MES_RealEmpty_';            % Output filename for Empty

% -----------------------------------------------------------

FileRecName = 'PCB';
Nbmic = 16;   % Number of microphones pairs (Total 32 mics)
nFile = 50;   % Number of acquisitions (Manual says 50)
ii0 = 1;      % Index of first test to load
FlagSave = 1; % 1 = Save data in mat format

%% ************************DON'T CHANGE NEXT PART*************************
%------------------------------
%% Initialize variables.
%------------------------------

%Acquisition sample parameters
ntot=8388608; % number of total bit sample
fe=12.5e6;    % sampling frequency of acquisition
fc=fe/2;      % sampling frequency of the clock
nbit = 24;    % number of bit for coding

%-------------------------------
% Preliminary calculations
%-------------------------------
nf=fe/fc; 
np=ntot/nf/64/2; 
np=np-1; 
np=floor(np/2)*2; 
fsample=fc/64/2;
tacq=np/fsample;
dt=1/fsample;

%%------------------------------
%Load I2S data in Time domain
%-------------------------------

%Initialization data matrices
yl=zeros(Nbmic,np);
yr=zeros(Nbmic,np);

disp(['Processing data from: ', FileDir]);
disp(['Saving to: ', SaveDir]);

for ii=ii0:1:ii0+nFile-1  % Loop on # of data acquisition
    tic;

    for mic=1:1:Nbmic % Loop over microphone pairs
        indexmic=mic;
        FName=strcat(FileRecName,num2str(indexmic),'Acq');
        
        % Construct file paths
        FileName = strcat(FileDir, FName, num2str(ii), '.csv');
        FileName_t = strcat(FileDir, FName, num2str(ii), '_time.csv');
        
        % Check if file exists before trying to read
        if isfile(FileName)
            %------------------------------
            % Data extraction
            %------------------------------       
            [r_val,l_val,t]=FCT_LoadData_I2SAcquisition_V3(FileName,FileName_t,np,nbit);
            yr(mic,:)=r_val(1:np); clear r_val
            yl(mic,:)=l_val(1:np); clear l_val    
            t=t(1:np);
        else
            warning(['File not found: ', FileName]);
        end
         
    end
    
    %------------------------------
    %Data organisation in a structure variable
    %------------------------------
    A.nbtot=ntot;
    A.f_acqu=fe;
    A.f_clk=fc;
    A.nbit=24;
    A.FileName=FileName;
    A.fsample=fsample;
    A.duration=tacq;
    A.Nb_real=nFile;
    A.ii_real=ii;
    A.Nb_mic=Nbmic;
    A.t=t;
    A.yr=yr;
    A.yl=yl;

    %------------------------------
    % Data saving
    %------------------------------
    if FlagSave==1
        % Check if Save Directory exists, if not create it
        if ~exist(SaveDir, 'dir')
           mkdir(SaveDir)
        end
        
        FileSave=strcat(SaveDir,FileSaveName,num2str(ii),'.mat');
        save(FileSave,'A');
    end
    sprintf('Antenna Acq n°%s, done in %f seconds',num2str(ii),toc)
    clear A

end

disp('Data loading complete!');