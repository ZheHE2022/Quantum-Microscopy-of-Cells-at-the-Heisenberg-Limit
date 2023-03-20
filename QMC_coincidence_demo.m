
%% This code implements an demo algorithm used in 
% "Quantum Microscopy of Cancer Cells at the Heisenberg Limit"

% Author: Zhe He
% Date: 01/12/2023

clear; clc; close all;
import matlab.io.*;

saveDataFolder = 'Coincidence/';

dataFolder = '';

filePrefix1 = "M1_USAF_bin2_group6_z1_12_X";
% centerX_left = (94+1)/2;
% centerY_left = (122+1)/2;
centerX_left = (98+1)/2;
centerY_left = (118+1)/2;
Xl = 0; Xr = 0; Yl = 0; Yr = 0;

center_name = ['_Xl',num2str(Xl), '_Xr',num2str(Xr), ...
    '_Yl',num2str(Yl), '_Yr',num2str(Yr)];


yideCode = true; % false: original software; true: LabVIEW

% make sure ceil and floor have difference
centerX_left = centerX_left + 0.01;
centerY_left = centerY_left + 0.01;


% plot during processing?
% plot_progress = true;
plot_progress = false;

%% two rounds
for filePrefix = [filePrefix1]
    
    %% Specify Files
    fullsavedataFolderPath = strcat(dataFolder, saveDataFolder);
    [status, msg, msgID] = mkdir(fullsavedataFolderPath);
    
    myFitFiles = dir(fullfile(dataFolder, strcat(filePrefix,'*.fits'))); %gets all txt files in struct
    for k = 1:length(myFitFiles)
        img_files{k}  = myFitFiles(k).name;
    end
    img_files = sort_nat(img_files,'ascend');
    
    spdc_photonsCurrent = readFile(strcat(dataFolder, img_files{1}), yideCode);
    arraySize = size(spdc_photonsCurrent);
    Nx = arraySize(1); Ny = arraySize(2);
    totalFrames = 0;
    
    %% define variables
    % Split the Cones Into Two Matrices of Equal Sizes Around the Middle
    % X is vertical (row); Y is horizontal (column)
    xHeight = min([Nx - ceil(centerX_left), floor(centerX_left) - 1]);
    xHeight = xHeight - 1;
    
    yWidth = min([Ny - ceil(centerY_left), floor(centerY_left) - 1]);
    yWidth = yWidth - 6;
    
    for total_data = 5
        finalCorMap_L = zeros(xHeight*2, yWidth);
        finalCorMap_R = finalCorMap_L;
        mean_leftSample = zeros(xHeight*2, yWidth);
        
        % adjust center based on 9x9 centers results
        centerX_left = centerX_left + Xl;
        centerY_left = centerY_left + Yl;
        
        batch_data = 5;
        for fileIndex0 = 0:batch_data:total_data-batch_data
            leftSample_list = [];
            aligned_Reference_list = [];
            
            disp(['Calculating Correlation Map for File: ', img_files{fileIndex0+1}]);
            for fileIndex = 1:batch_data
                
                spdc_photons = readFile(strcat(dataFolder, img_files{fileIndex0+fileIndex}), yideCode);
                
                if plot_progress && (fileIndex == 1)
                    % Plot the Frame Average
                    frame_ave = mean(spdc_photonsCurrent,3); % Result is a y,x dimensional array
                    figure(1);
                    set(gcf, 'Position', [0, 0, 500, 400])
                    imagesc(frame_ave, [min(frame_ave(:)), max(frame_ave(:))])
                    axis equal tight
                    colormap hot
                    colorbar
                    title(strcat("Frame average of file #", num2str(fileIndex)));
                    drawnow
                end
                
                %% Get the Correlation Map from Photons
                
                
                leftSample = spdc_photons(ceil(centerX_left) - xHeight:floor(centerX_left) + xHeight, ceil(centerY_left)-yWidth:floor(centerY_left), :);
                leftReference = spdc_photons(ceil(centerX_left) - xHeight:floor(centerX_left) + xHeight, ceil(centerY_left):floor(centerY_left) + yWidth, :);
                
                
                
                % Align the Right Values with the Left
                aligned_Reference = flip(flip(leftReference, 2), 1);
           
                % Plot the Left and Right Split Side by Side    
                mean_Reference = mean(aligned_Reference,3);
                
                if plot_progress && (fileIndex == 1)
                    figure(2)
                    set(gcf, 'Position', [500, 0, 900, 400])
                    
                    subplot(1,3,1)
                    imagesc(mean_leftSample)
                    axis equal tight
                    colormap hot
                    colorbar
                    title('LeftSample');
                    
                    subplot(1,3,2)
                    imagesc(mean_Reference)
                    axis equal tight
                    colormap hot
                    colorbar
                    title('Reference');
                    
                    
                    subplot(1,3,3)
                    imagesc(imfuse(mean_Reference,mean_leftSample,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]))
                    axis equal tight
                    title("Left Fuse")
                    
                    
                    
                    
                    exportgraphics(gcf, strcat(fullsavedataFolderPath, filePrefix, center_name, '_Align.tif'), 'Resolution', '300');
                    
                    %         drawnow
                end
                
                %% accumulate the matrices
                leftSample_list = cat(3, leftSample_list, leftSample);
                aligned_Reference_list = cat(3, aligned_Reference_list, aligned_Reference);
                                 
            end
            
            %% calculate the coincidence
            corMapCurrent_Left = corrcoin(leftSample_list, aligned_Reference_list);
            
            %% Add to Running Weighted Average of Correlation Map
            
            mean_leftSample = mean_leftSample + mean(leftSample_list, 3);
            finalCorMap_L = finalCorMap_L + (corMapCurrent_Left);
        end
        
        % Average Out All the Frames
        totalFrames = total_data/batch_data;
        mean_leftSample = mean_leftSample./totalFrames;
        finalCorMap_L = finalCorMap_L./totalFrames;
        
        %% Plot Correlation Map
        
        figure(4);
        set(gcf, 'Position', [500, 400, 600, 300])
        imagesc((finalCorMap_L))
        axis equal tight
        colormap hot
        xlabel('K _X_1 + K _X_2')
        ylabel('K _Y_1 + K _Y_2')
        colorbar
        caxis([0, 80])
        title('QMC')

        exportgraphics(gcf, strcat(fullsavedataFolderPath, filePrefix, center_name, '_CoincidenceMap_', num2str(total_data), '.tif')...
            , 'Resolution', '300');
        
        save(strcat(fullsavedataFolderPath, filePrefix, center_name, '_coincidencedata_X', num2str((centerX_left-0.01)*2-1)...
            , '_Y', num2str((centerY_left-0.01)*2-1), '_', num2str(total_data), '.mat'),...
            'finalCorMap_L','mean_leftSample');
        

    end
end


%% Define Function
function spdc_photons0 = readFile(spdc_Photons_File, yideCode)
import matlab.io.*;
file_Pointer_spdc = fits.openFile(spdc_Photons_File);
if yideCode
    fits.movAbsHDU(file_Pointer_spdc, 2);
end
spdc_photons0 = single(fits.readImg(file_Pointer_spdc));
fits.closeFile(file_Pointer_spdc);
end



function [cs,index] = sort_nat(c,mode)
if nargin < 2
    mode = 'ascend';
end
% Make sure mode is either 'ascend' or 'descend'.
modes = strcmpi(mode,{'ascend','descend'});
is_descend = modes(2);
if ~any(modes)
    error('sort_nat:sortDirection',...
        'sorting direction must be ''ascend'' or ''descend''.')
end
% Replace runs of digits with '0'.
c2 = regexprep(c,'\d+','0');
% Compute char version of c2 and locations of zeros.
s1 = char(c2);
z = s1 == '0';
% Extract the runs of digits and their start and end indices.
[digruns,first,last] = regexp(c,'\d+','match','start','end');
% Create matrix of numerical values of runs of digits and a matrix of the
% number of digits in each run.
num_str = length(c);
max_len = size(s1,2);
num_val = NaN(num_str,max_len);
num_dig = NaN(num_str,max_len);
for i = 1:num_str
    num_val(i,z(i,:)) = sscanf(sprintf('%s ',digruns{i}{:}),'%f');
    num_dig(i,z(i,:)) = last{i} - first{i} + 1;
end
% Find columns that have at least one non-NaN.  Make sure activecols is a
% 1-by-n vector even if n = 0.
activecols = reshape(find(~all(isnan(num_val))),1,[]);
n = length(activecols);
% Compute which columns in the composite matrix get the numbers.
numcols = activecols + (1:2:2*n);
% Compute which columns in the composite matrix get the number of digits.
ndigcols = numcols + 1;
% Compute which columns in the composite matrix get chars.
charcols = true(1,max_len + 2*n);
charcols(numcols) = false;
charcols(ndigcols) = false;
% Create and fill composite matrix, comp.
comp = zeros(num_str,max_len + 2*n);
comp(:,charcols) = double(s1);
comp(:,numcols) = num_val(:,activecols);
comp(:,ndigcols) = num_dig(:,activecols);
% Sort rows of composite matrix and use index to sort c in ascending or
% descending order, depending on mode.
[~,index] = sortrows(comp);
if is_descend
    index = index(end:-1:1);
end
index = reshape(index,size(c));
cs = c(index);
end

function corrcoin = corrcoin(Left, Right)
meanLeft = mean(Left, 3);
meanRight = mean(Right, 3);
corrcoin = mean((Left - meanLeft) .* (Right - meanRight), 3);
end
