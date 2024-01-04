function [artifacts] = detect_artifacts_hjorth(data, Fs, window_length_sec)
%EEG_DETECT_TIME_DOMAIN_ARTIFACTS  Detect artifacts in the time domain by iteratively removing data above a given z-score criterion
%
%   Usage:
%   Direct input:
%   artifacts = detect_artifacts(data, Fs, crit_units, hf_crit, hf_pass, bb_crit, bb_pass, smooth_duration, ...
%                                            verbose, histogram_plot, return_filts_only, hpFilt_high, hpFilt_broad, detrend_filt)
%
%   Input:
%   data: 1 x <number of samples> vector - time series data-- required
%   Fs: double - sampling frequency in Hz  -- required
%   crit_units: string 'std' to use iterative crit_units (default) or a strict threshold on 'MAD',defined as K*MEDIAN(ABS(A-MEDIAN(A)))
%   hf_crit: double - high frequency criterion - number of stds/MAD above the mean to remove (default: 3.5)
%   hf_pass: double - high frequency pass band - frequency for high pass filter in Hz (default: 25 Hz)
%   bb_crit: double - broadband criterion - number of stds/MAD above the mean to remove (default: 3.5)
%   bb_pass: double - broadband pass band - frequency for high pass filter in Hz (default: .1 Hz)
%   smooth_duration: double - time (in seconds) to smooth the time series (default: 2 seconds)
%   verbose: logical - verbose output (default: false)
%   histogram_plot: logical - plot histograms for debugging (default: false)
%   return_filts_only: logical - return 3 digitalFilter objects and nothing else [1x3] (order is hpFilt_high, hpFilt_broad, detrend_filt)
%                                for use in this function (default: false)
%   hpFilt_high: digitalFilter - Includes parameters to use for the high frequency high pass filter
%   hpFilt_broad: digitalFilter - Includes parameters to use for the broadband high pass filter
%   detrend_filter: digitalFilter - Includes parameters to use for the detrending high pass filter
%
%   Output:
%   artifacts: 1xT logical of times flagged as artifacts (logical OR of hf and bb artifacts)
%
%   Copyright 2020 Michael J. Prerau, Ph.D. - http://www.sleepEEG.org
%   This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
%   (http://creativecommons.org/licenses/by-nc-sa/4.0/)
%
%   Last modified 01/22/2021 by Mike & Tom, Alex
%% ********************************************************************

%% Input parsing
%Force column vector for uniformity
if ~iscolumn(data)
    data=data(:);
end

%Force to double
if ~isa(data,'double')
    data = double(data);
end

if nargin<2
    error('Data vector and sampling rate required');
end

if nargin<3
    window_length_sec = 4;
end 

%% Get bad indicies
%Get bad indices
bad_inds = (isnan(data) | isinf(data) | find_flat(data))';

%Interpolate big gaps in data
tt = 1:length(data);
data_fixed = interp1([0, tt(~bad_inds), length(data)+1], [0; data(~bad_inds); 0], tt)';

%% Compute hjorth parameters for overlapping windows  
% Using linear filters for faster computation
[ activity, mobility, complexity ] = hjorth(data_fixed, round(window_length_sec * Fs));

% % Manual for loop to check implementation
% % Set up segments 
% window_params(1) = window_length_sec; window_params(2) = 0.05;
% window_length_samples = round(window_params(1) * Fs);
% window_step_samples = round(window_params(2) * Fs);
% begin_idx = (1:window_step_samples:length(data_fixed)-window_length_samples+1)';
% window_center_secs = (begin_idx+window_length_samples/2-1)/Fs;
% 
% activity = zeros(size(begin_idx));
% mobility = zeros(size(begin_idx));
% complexity = zeros(size(begin_idx));
% 
% for ii = 1:length(begin_idx)
%     y = data_fixed(begin_idx(ii):begin_idx(ii)+window_length_samples-1);
%     dy = diff([0;y],1);
%     ddy = diff([0;dy],1);
%     
%     var_y = var(y);
%     var_dy = var(dy);
%     var_ddy = var(ddy);
%     
%     activity(ii) = var_y;
%     mobility(ii) = sqrt(var_dy/var_y);
%     complexity(ii) = sqrt(var_ddy/var_dy)/mobility(ii);
% end

%% DEBUGGING PLOT
figure;
ax = figdesign(6,1, 'merge', {2:3},'type','usletter','orient','landscape');

axes(ax(1))
hypnoplot(t, stages.stages)

axes(ax(2))
imagesc(stimes, sfreqs, nanpow2db(spect));
climscale;
axis xy
xlabel('Time (s)');
ylabel('Frequency (Hz)');
cbar = topcolorbar;
cbar.Position(2) = cbar.Position(2) - .05; %increase the number to move it further down
colormap(ax(2), jet);
axis tight
title('Full-night Spectrogram'); 
set(gca,'FontSize',8)

hline(17.25);
hline(7.75);

axes(ax(3))
plot(t, activity)
xlabel('Activity')

axes(ax(4))
plot(t, mobility)
xlabel('Mobility')

axes(ax(5))
plot(t, complexity)
xlabel('Complexity')

linkaxes(ax, 'x');

%% Join artifacts from different frequency bands (not yet updated)
% artifacts = hf_artifacts | bb_artifacts;
% % sanity check before outputting
% assert(length(artifacts) == length(data), 'Data vector length is inconsistent. Please check.')

%Find all the flat areas in the data
function binds = find_flat(data, min_size)
if nargin<2
    min_size = 100;
end

%Get consecutive values equal values
[clen, cind] = getchunks(data);

%Return indices
if isempty(clen)
    inds = [];
else
    size_inds = clen>=min_size;
    clen = clen(size_inds);
    cind = cind(size_inds);
    
    flat_inds = cell(1,length(clen));
    
    for ii = 1:length(clen)
        flat_inds{ii} = cind(ii):(cind(ii)+(clen(ii)-1));
    end
    
    inds = cat(2,flat_inds{:});
end

binds = false(size(data));
binds(inds) = true;
