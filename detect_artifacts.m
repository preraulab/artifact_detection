function artifacts = detect_artifacts(data, Fs, varargin)
%DETECT_ARTIFACTS  Detect and remove artifacts in time series data
%
%   Usage:
%       artifacts = detect_artifacts(data, Fs, varargin)
%   
%   NOTE: Default arguments will be the same as the original
%   detect_artifacts()
%
%   Input:
%       data: <number of samples> x 1 vector - time series data -- required
%       Fs: double - sampling frequency in Hz -- required
%
%   Optional Input Arguments:
%       'isexcluded': logical vector - indicates excluded data points that will be set as artifact (default: [])
%       'zscore_method': char - method for z-score calculation ('standard':mean and std or 'robust': median and MAD) (default: 'standard')
%       'hf_crit': numeric - high-frequency artifact detection criterion (default: 4.5)
%       'hf_pass': numeric - high-frequency artifact passband frequency in Hz (default: 35)
%       'bb_crit': numeric - broadband artifact detection criterion (default: 4.5)
%       'bb_pass': numeric - broadband artifact passband frequency in Hz (default: 0.1)
%       'hf_detrend': logical - apply detrending to high-frequency artifacts (default: true)
%       'bb_detrend': logical - apply detrending to broadband artifacts (default: true)
%       'slope_test': logical - enable slope test for artifacts (default: false)
%       'slope_crit': numeric - slope threshold for the slope test (default: -0.5)
%       'smooth_duration': numeric - duration for data smoothing in seconds (default: 2)
%       'detrend_duration': numeric - duration for detrending in seconds (default: 300)
%       'buffer_duration': numeric - duration for adding buffer around detected artifacts in seconds (default: 0)
%       'verbose': logical - enable verbose mode (default: false)
%       'diagnostic_plot': logical - enable diagnostic plots (default: false)
%       'histogram_plot': logical - enable histogram plot during artifact detection (default: false)
%       'return_filts_only': logical - return filter parameters only (default: false)
%       'hpFilt_high': [] - high-pass filter design (default: [])
%       'hpFilt_broad': [] - broadband high-pass filter design (default: [])
%
%   Output:
%       artifacts: logical vector - indicates artifact time points
%
%   Example:
%   To detect and remove artifacts in EEG data:
%       Fs = 1000; % Sampling Frequency
%       data = randn(10000, 1); % Example EEG data
%
%       artifact = detect_artifacts_dev(data,Fs);
%
%       %Old version
%       artifacts_old = detect_artifacts_dev(data, Fs, 'zscore_method','standard','hf_crit', 4.5,'bb_crit', 4.5,'slope_test',false);
%
%   Copyright 2024 Michael J. Prerau Laboratory. - http://www.sleepEEG.org
%
%   See also: multitaper_spectrogram_mex, get_chunks, consecutive_runs
%% ********************************************************************

%% Input parsing
%Force column vector for uniformity
if ~iscolumn(data)
    data=data(:);
end

%Force to double
if ~isa(data, 'double')
    data = double(data);
end

p = inputParser;

addOptional(p, 'isexcluded', [], @(x) validateattributes(x, {'logical', 'vector'},{}));
addOptional(p, 'zscore_method', 'robust',  @(x) validateattributes(x, {'char'},{}));

addOptional(p, 'hf_crit', 5.5, @(x) validateattributes(x,{'numeric'},{'real', 'finite', 'nonnan'}));
addOptional(p, 'hf_pass', 35, @(x) validateattributes(x,{'numeric'},{'real', 'finite', 'nonnan'}));
addOptional(p, 'bb_crit', 5.5, @(x) validateattributes(x,{'numeric'},{'real', 'finite', 'nonnan'}));
addOptional(p, 'bb_pass', .1, @(x) validateattributes(x,{'numeric'},{'real', 'finite', 'nonnan'}));
addOptional(p, 'hf_detrend', true, @(x) validateattributes(x,{'logical'},{'real','nonempty', 'nonnan'}));
addOptional(p, 'bb_detrend', true, @(x) validateattributes(x,{'logical'},{'real','nonempty', 'nonnan'}));

addOptional(p, 'slope_test', true, @(x) validateattributes(x,{'logical'},{'real','nonempty', 'nonnan'}));
addOptional(p, 'slope_crit', -0.5, @(x) validateattributes(x,{'numeric'},{'real', 'finite', 'nonnan'}));

addOptional(p, 'smooth_duration', 2, @(x) validateattributes(x,{'numeric'},{'real', 'finite', 'nonnan'}));
addOptional(p, 'detrend_duration', 5*60, @(x) validateattributes(x,{'numeric'},{'real', 'finite', 'nonnan'}));
addOptional(p, 'buffer_duration', 0, @(x) validateattributes(x,{'numeric'},{'real', 'finite', 'nonnan'}));

addOptional(p, 'verbose', false, @(x) validateattributes(x,{'logical'},{'real','nonempty', 'nonnan'}));
addOptional(p, 'diagnostic_plot', false, @(x) validateattributes(x,{'logical'},{'real','nonempty', 'nonnan'}));
addOptional(p, 'histogram_plot', false, @(x) validateattributes(x,{'logical'},{'real','nonempty', 'nonnan'}));
addOptional(p, 'return_filts_only', false, @(x) validateattributes(x,{'logical'},{'real','nonempty', 'nonnan'}));

addOptional(p, 'hpFilt_high', []);
addOptional(p, 'hpFilt_broad', []);

parse(p,varargin{:});
parser_results = struct2cell(p.Results); %#ok<NASGU>
field_names = fieldnames(p.Results);

eval(['[', sprintf('%s ', field_names{:}), '] = deal(parser_results{:});']);

validatestring(zscore_method,{'robust','standard'});

% Verify that the isexcluded boolean is valid
if isempty(isexcluded)
    isexcluded = false(size(data));
else
    if ~iscolumn(isexcluded)
        isexcluded=isexcluded(:);
    end
    assert(length(isexcluded) == length(data),'isexcluded must be the same length as data');
end

% Create filters if none are provided
if isempty(hpFilt_high) %#ok<*NODEF>
    hpFilt_high = designfilt('highpassiir','FilterOrder',4, ...
        'PassbandFrequency',hf_pass,'PassbandRipple',0.2, ...
        'SampleRate',Fs);
end

if isempty(hpFilt_broad)
    hpFilt_broad = designfilt('highpassiir','FilterOrder',4, ...
        'PassbandFrequency',bb_pass,'PassbandRipple',0.2, ...
        'SampleRate',Fs);
end

% If desired, return filter parameters only
if return_filts_only
    artifacts = [hpFilt_high, hpFilt_broad];
    return
end

if verbose
    t_start = tic;
    disp('Performing artifact detection:')
    disp(['     Z-score method: ' zscore_method])
    disp(' ')
    disp(['     High Frequency Criterion: ' num2str(hf_crit)]);
    disp(['     High Frequency Passband: ' num2str(hf_pass) ' Hz']);
    disp(['     High Frequency Detrend: '  char(string(hf_detrend))]);
    disp(' ')
    disp(['     Broadband Criterion: ' num2str(bb_crit) ]);
    disp(['     Broadband Passband: ' num2str(bb_pass) ' Hz']);
    disp(['     Broadband Detrend: '  char(string(hf_detrend))]);
    disp(' ')

    if slope_test
        disp(['     Slope Test: slope > ' num2str(slope_crit)]);
    else
        disp('     Slope Test: False');
    end
    disp('    ');
end

%% Get bad indicies

if slope_test
    [spect, stimes, sfreqs] = multitaper_spectrogram_mex(data, Fs, [1 min(55,Fs/2)],[10,19],[10 5],[],[],[],false,false);
    B = zeros(length(stimes),2);
    for ii = 1:length(stimes)
        B(ii,:) = polyfit(log(sfreqs)',log(spect(:,ii)),1);
    end
    t = (0:length(data)-1)/Fs;
    bad_slope = interp1(stimes,double(B(:,1)>slope_crit)',t,'nearest')';
    bad_slope(isnan(bad_slope)) = 1;
else
    bad_slope = zeros(size(data));
end

%Get bad indices
[~, ~, ~, is_flat] = get_chunks(data, Fs);
bad_inds = isnan(data) | isinf(data) | is_flat | bad_slope | isexcluded;
bad_inds = find_outlier_noise(data, bad_inds);

%Interpolate big gaps in data
t = 1:length(data);
good_data = data(~bad_inds);
data_fixed = interp1([0, t(~bad_inds), length(data)+1], [good_data(1); good_data; good_data(end)], t)';

%Trim the data to make it more consistent with cutting vs excluding
valid_range_inds = find(~bad_inds,1,'first'):find(~bad_inds,1,'last');
data_fixed = data_fixed(valid_range_inds);
bad_inds = bad_inds(valid_range_inds);

%% Get high frequency artifacts
hf_artifacts = compute_artifacts(hpFilt_high, detrend_duration, hf_crit, data_fixed, smooth_duration, Fs, bad_inds, verbose,...
    'high frequency', histogram_plot, zscore_method, hf_detrend, diagnostic_plot);

%% Get broad band frequency artifacts
bb_artifacts = compute_artifacts(hpFilt_broad, detrend_duration, bb_crit, data_fixed, smooth_duration, Fs, bad_inds, verbose,...
    'broadband frequency', histogram_plot, zscore_method, bb_detrend, diagnostic_plot);

%% Join artifacts from different frequency bands and put back to original size
artifacts = true(size(data));
artifacts(valid_range_inds) = hf_artifacts | bb_artifacts | bad_inds;

%% Add buffer on both sides of detected artifacts
if buffer_duration > 0
    [cons, inds] = consecutive_runs(artifacts);
    for ii = 1:length(cons)
        buffer_start_idx = ceil(max(1, inds{ii}(1)-buffer_duration*Fs));
        buffer_end_idx = ceil(min(length(artifacts), inds{ii}(end)+buffer_duration*Fs));
        artifacts(buffer_start_idx:buffer_end_idx) = true;
    end
end

%% Sanity check before outputting
assert(length(artifacts) == length(data), 'Data vector length is inconsistent. Please check.')

if verbose
    disp(['Total time: ' num2str(toc(t_start)) ' seconds']);
    disp(' ');
end

%Find all time points that have outlier high or low noise
function bad_inds = find_outlier_noise(data, bad_inds)
data_mean = mean(data(~bad_inds));
data_std = std(data(~bad_inds));

outlier_scalar = 10;
low_thresh = data_mean - outlier_scalar * data_std;
high_thresh = data_mean + outlier_scalar * data_std;

inds = data <= low_thresh | data >= high_thresh;
bad_inds(inds) = true;


function detected_artifacts = compute_artifacts(filter_coeff, detrend_duration, crit, data_fixed,...
    smooth_duration, Fs, bad_inds, verbose, verbosestring, histogram_plot, zscore_method, detrend_on, diagnostic_plot)
%% Get artifacts for a particular frequency band

if diagnostic_plot
    figure
    if detrend_on
        ax = figdesign(5,1);
    else
        ax = figdesign(4,1);
    end
    linkaxes(ax,'x');
    t = (0:length(data_fixed)-1)/Fs;
end

%Perform a high pass filter
y_signal = filter(filter_coeff, data_fixed);

if diagnostic_plot
    axes(ax(1))
    plot(t,y_signal)
    title('Filtered Signal')
end

%Look at the data envelope
y_signal = abs(hilbert(y_signal));

if diagnostic_plot
    axes(ax(2))
    plot(t,y_signal)
    title('Hilbert Magnitude')
end

%Smooth data
y_signal = movmean(y_signal, smooth_duration*Fs);

if diagnostic_plot
    axes(ax(3))
    plot(t,y_signal)
    title('Smoothed Magnitude')
end

% We should smooth then take log
y_signal = log(y_signal);

if diagnostic_plot
    axes(ax(4))
    plot(t,y_signal)
    title('Log Smoothed Magnitude')
end

% Detrend data
if detrend_on
    %Moving median is the most efficient way to detrend and does not leave
    %a transient
    y_signal = y_signal - movmedian(y_signal, Fs*detrend_duration);
    % y_detrend = spline_detrend(y_log, Fs, [], 60);
    % y_detrend = filter(detrend_filt, y_log);

    if diagnostic_plot
        axes(ax(5))
        plot(t,y_signal)
        title('Detrended Log Smoothed Magnitude');
    end
end

% create an indexing array marking bad timepoints before z-scoring
detected_artifacts = bad_inds;

%Take z-score
ysig = y_signal(~detected_artifacts);
if strcmpi(zscore_method,'standard')
    ymid = mean(ysig);
    ystd = std(ysig);
else
    ymid = median(ysig);
    ystd = mad(ysig);
end

y_signal = (y_signal - ymid)/ystd;

if verbose
    num_iters = 1;
    disp(['     Running ', verbosestring, ' detection...']);
end

if histogram_plot
    fh = figure;
    set(gca,'nextplot','replacechildren');
    histogram(y_signal,100);
    title(['          Iteration: ' num2str(num_iters)]);
    drawnow;
    ah = gca;
end


%Keep removing until all values under criterion
over_crit = abs(y_signal)>crit & ~detected_artifacts;

%Loop until nothing over criterion
count = 1;
while any(over_crit)

    %Update the detected artifact time points
    detected_artifacts(over_crit & ~bad_inds) = true;

    %Compute modified z-score
    ysig = y_signal(~detected_artifacts);
    if strcmpi(zscore_method,'standard')
        ymid = mean(ysig);
        ystd = std(ysig);
    else
        ymid = median(ysig);
        ystd = mad(ysig);
    end

    y_signal = (y_signal - ymid)/ystd;

    %Find new criterion
    over_crit = abs(y_signal)>crit & ~detected_artifacts;

    if histogram_plot
        axes(ah); %#ok<LAXES>
        histogram(y_signal(~detected_artifacts), 100);
        title(['          Outliers Removed: iteration ', num2str(count)]);
        drawnow;
        pause(0.1);
    end
    count = count + 1;
end

if verbose
    disp(['     Ran ' num2str(num_iters) ' iterations']);
end

if histogram_plot
    close(fh);
end
