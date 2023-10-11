%% Three-Point Bending Analysis

clear all
close all
clc

%{
Updated 13 September 2023
Author: Macy Mora-Antoinette

MATLAB R2019a  

This matlab code is intended to help users compute the mechanical
properties from force-displacement data collected typically by three-point
bending of mouse femurs. 

It is separate into two parts:
PART 1: A simple for loop to iterate through all the raw data files
outputted by the MTS (mechanical testing system) located in Biomechanics
Labs (I believe the machine is owned by the van der Meulen group). It will
save the results in a csv.
PART 2: This function is called in PART 1 and is the main function that
computes the mechanical properties. It will save the plots of the
force-displacement data, which is good to have for the record. 

Key things to keep in mind:
1. This code is not automated. It requires user input for every data file.
The user will select the datapoints for stiffness, yeild force, etc. 
2. You must change the force to voltage equation in line 64. Typically,
a new calibration curve is created just before a round of testing, and you
try to use the same curve for the entire experiment by testing all your 
bones within a few days.
3. The columns for displacement and voltage can change depending on the
testing regime you pick on the MTS, so double check they're correct.
4. The FileName is the name you type in the MTS software when testing a new
specimen. This code requires a specific naming format. 
Ex C21-R-F-RF (mouseline-eartag-sex-rightfemur) all hyphens, NO SPACES
You can change this is you change the splitting mechanisms in lines 57 and 196

                                --- PART 1 ---
%}
%%
% --------------  IMPORT DATA -------------- %
columnDisp = 1; %raw data collects displacement in column 1
columnVolt = 4; %raw data collects voltage output in column 4

%save data directory - make sure ONLY DATA FOLDERS are in this directory
disp('Select Directory that has only Data folders')
DatDir = uigetdir(); 
ls = dir(DatDir); %get list of contents
dat = {ls.name}; %all folders listed

Results = {};
%Iterate over all datafolders
for i=3:length(dat) %first two elements are extra, skip
    dat_folder = dat{i}; 
    FileName = split(dat_folder, " "); FileName = FileName{1}; %name of sample
    resultsDir = [DatDir, '\', dat_folder];
    data_file = dir(resultsDir); data_file = {data_file.name};
    data = readmatrix([resultsDir, '\', data_file{3}]); %import data

    displacement = -1.*data(:,columnDisp);%column with mm data, changes orientation of graph to familiar one
    voltage = -1.*data(:,columnVolt);%column with volt data, changes orientation of graph to familiar one
    force = 72.572.*voltage + 0.0552;%converts from Volts to Newtons, take from load cell calibration curve
    
    %Table of mechanical properties
    outputs = mechanical_properties(displacement, force, FileName);
    
    %make table with all outputs as mechanical properties
    Results = [Results;outputs];
end

%create table to store results
T = cell2table(Results); %empty table
T.Properties.VariableNames = {'Physical_Tag', 'Alt_ID', 'Sex', 'Femur', 'max_force', 'ultimate_energy',...
    'failure_energy', 'fracture_force', 'fracture_disp',...
    'stiffness', 'yield_force', 'yield_disp', 'PYD', 'yield_energy'};
%create csv files
writetable(T, 'Three_Point_Bending_Data.csv');

%%
%{
                                --- PART 2 ---
%}

function outputs = mechanical_properties(displacement, force, FileName)
    
% Description: inputs are force and displacement. The output is an
% array of all of the mechanical properties computed from the
% force-displacement curve
     
    % --------------  PREPROCESS DATA -------------- %

    
    %During our experiment, we did not preload the bone, so the first several
    %seconds we only see noise. Identify the expected start true loading on the bone
    %An alternative is to preload the bone and avoid this step, but be consistent
    %within a single experiment
    x = displacement;
    y = force;
    plot(displacement, force)
    xlabel('Deflection (mm)')
    ylabel('Load (N)')
    title('Remove nonloading region (Identify start and break)');
    [xinput,~] = ginput(2); %user selects two point on the plot
    %Select data points closest to the points selected
    [~, data_start] = min(abs(x-xinput(1)));
    [~, data_end] = min(abs(x-xinput(2)));
    x = x(data_start:data_end);
    y = y(data_start:data_end);
    x = x - x(1);
    
    %remove noise from data to help compute parameters
    smoothy = smooth(y);
    smoothx = smooth(x);
    
    % --------------  COMPUTE MECHANICAL PROPERTIES -------------- %
    % Mechanical Properties
    % 1. Max or Ultimate Force/Load and Energy
    % 2. Work to Fracture/Energy to Failure
    % 3. Fracture or Failure Force/Load and Displacement
    % 4. Stiffness
    % 5. Yield Force/Load
    % 6. Post-yield Displacement and Yield Displacement
    % 7. Yield Energy 

    % 1. max/ultimate force - max of entire curve
    [~, maxF_idx] = max(smoothy);
    max_force = y(maxF_idx);
    ultimate_energy = trapz(x(1:maxF_idx), y(1:maxF_idx));

    % 2. work to fracture (Energy to failure) - area under the curve
    failure_energy = trapz(x, y);

    % 3. Fracture/breaking force or load - end of curve
    fracture_force = y(end);
    fracture_disp = x(end);

    % 4. stiffness - slope of linear region
    plot(smoothx, smoothy)
    xlabel('Deflection (mm)')
    ylabel('Load (N)')
    title('Click on 2 Points (Left to Right) in the Linear Region for the Slope');
    [xinput,yinput] = ginput(2); %user selects two point on the plot
    while xinput(1) > xinput(2)  %makes sure second click is up and to the right of the first click
        title('Invalid Selection: Select 2 Points (Left to Right) in the Linear Region')
        [xinput,yinput] = ginput(2); %user selects two point on the plot
    end
    %Compute stiffness based on the points closest to the points selected
    [~, Start] = min(abs(x-xinput(1)));
    [~, End] = min(abs(x-xinput(2)));
    deltay = smoothy(End) - smoothy(Start);
    deltax = smoothx(End) - smoothx(Start);
    stiffness = deltay/deltax;
    stiffness_eqn = stiffness.*(smoothx - smoothx(Start)) + smoothy(Start);

    % 5. Yield Strength
    %find the beginning of the elatic region by finding the closest
    %intersection to the stiffness (linear) equation
    dist = abs(y - stiffness_eqn); %absolute differences
    intersect = find(dist < 0.1); %consider 0.1 as the minimum distance to be considered an intersection
    intercept = intersect(1); %first intersection is the elastic region intercept
    
    %find yield strength by finding the closest
    %intersection to the yield (stiffness with 10% reduced slope) equation
    yield_stiffness = 0.9*stiffness; % 10% rule to find yield
    yield_eqn = yield_stiffness.*(smoothx - smoothx(intercept)) + smoothy(intercept);
    
    figure;
    plot(x,y,'r', smoothx, yield_eqn, 'k')
    xlabel('Deflection (mm)')
    ylabel('Load (N)')
    title('Select the Intersection for Yield Force')
    [x_input,y_input] = ginput(1);
    [~, yield_idx] = min(abs(x-x_input(1)));
    yield_force = y(yield_idx);
    yield_disp = x(yield_idx);

    % 6. Post-yield Displacement and Yield Displacement
    PYD = fracture_disp - yield_disp;

    % 7. Yield Energy (Energy to Yield) - Area under yield curve
    yield_energy = trapz(x(1:yield_idx), y(1:yield_idx));

    %Plot final results
    figure;
    plot(x,y,'r', smoothx, stiffness_eqn,'b', x(maxF_idx), max_force,'o', smoothx, yield_eqn, 'k', yield_disp, yield_force, 'o')
    xlabel('Deflection (mm)')
    ylabel('Load (N)')
    title('Stiffness, Yield Force, and Ultimate Force')
    saveas(gcf, FileName)
    close all
    
    %output arguments
    %Name encodes mouse information
    mouse = split(FileName, "-"); 
    %mechanical properties
    mech_prp = {max_force, ultimate_energy, failure_energy, fracture_force,...
        fracture_disp, stiffness, yield_force, yield_disp, PYD, yield_energy};
    
    outputs = [mouse' mech_prp];
end