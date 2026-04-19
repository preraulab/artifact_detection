function artifacts = detect_artifacts(data, Fs, varargin)
%DETECT_ARTIFACTS  Detect and remove artifacts in time series data
%
%   Usage:
%       artifacts = detect_artifacts(data, Fs, varargin)
%
%   Input:
%       data: <number of samples> x 1 vector - time series data -- required
%       Fs: double - sampling frequency in Hz -- required
%
%   Optional Input Arguments:
%       'isexcluded': logical vector - indicates excluded data points that will be set as artifact (default: logical([]))
%       'exclude_mode': string - 'data' or 'artifact', determines if isexcluded time are treated as data
%                       that are not factored in for the theshold detection or artifacts (default: 'data')
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
p = inputParser;

addOptional(p, 'isexcluded', logical([]), @(x) validateattributes(x,{'logical'},{'real','finite','2d'}));
addOptional(p, 'exclude_mode', 'data', @(x) any(validatestring(x, {'artifact', 'data'})));
addOptional(p, 'zscore_method', 'robust', @(x) any(validatestring(x, {'robust', 'standard'})));

addOptional(p, 'hf_crit', 5.5, @(x) validateattributes(x,{'numeric'},{'real','finite','scalar'}));
addOptional(p, 'hf_pass', 35, @(x) validateattributes(x,{'numeric'},{'real','finite','scalar'}));
addOptional(p, 'bb_crit', 5.5, @(x) validateattributes(x,{'numeric'},{'real','finite','scalar'}));
addOptional(p, 'bb_pass', .1, @(x) validateattributes(x,{'numeric'},{'real','finite','scalar'}));
addOptional(p, 'hf_detrend', true, @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addOptional(p, 'bb_detrend', true, @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

addOptional(p, 'slope_test', true, @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addOptional(p, 'slope_crit', -0.5, @(x) validateattributes(x,{'numeric'},{'real','finite','scalar'}));

addOptional(p, 'smooth_duration', 2, @(x) validateattributes(x,{'numeric'},{'real','finite','scalar'}));
addOptional(p, 'detrend_duration', 5*60, @(x) validateattributes(x,{'numeric'},{'real','finite','scalar'}));
addOptional(p, 'buffer_duration', 0, @(x) validateattributes(x,{'numeric'},{'real','finite','scalar'}));

addOptional(p, 'verbose', false, @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addOptional(p, 'diagnostic_plot', false, @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addOptional(p, 'histogram_plot', false, @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addOptional(p, 'return_filts_only', false, @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

addOptional(p, 'hpFilt_high', []);
addOptional(p, 'hpFilt_broad', []);

parse(p,varargin{:});

opts = p.Results;

%Force column vector for uniformity
if isrow(data)
    data = data(:);
end

% Verify that the isexcluded boolean is valid
if isempty(opts.isexcluded)
    opts.isexcluded = false(size(data));
else
    if isrow(opts.isexcluded)
        opts.isexcluded=opts.isexcluded(:);
    end
    assert(length(opts.isexcluded) == length(data),'isexcluded must be the same length as data');
end

% Create filters if none are provided
if isempty(opts.hpFilt_high) %#ok<*NODEF>
    opts.hpFilt_high = designfilt('highpassiir','FilterOrder',4, ...
        'PassbandFrequency',opts.hf_pass,'PassbandRipple',0.2, ...
        'SampleRate',Fs);
end

if isempty(opts.hpFilt_broad)
    opts.hpFilt_broad = designfilt('highpassiir','FilterOrder',4, ...
        'PassbandFrequency',opts.bb_pass,'PassbandRipple',0.2, ...
        'SampleRate',Fs);
end

% If desired, return filter parameters only
if opts.return_filts_only
    artifacts = [opts.hpFilt_high, opts.hpFilt_broad];
    return
end

if opts.verbose
    t_start = tic;
    disp('Performing artifact detection:')
    disp(['     Z-score method: ' opts.zscore_method])
    disp(' ')
    disp(['     High Frequency Criterion: ' num2str(opts.hf_crit)]);
    disp(['     High Frequency Passband: ' num2str(opts.hf_pass) ' Hz']);
    disp(['     High Frequency Detrend: '  char(string(opts.hf_detrend))]);
    disp(' ')
    disp(['     Broadband Criterion: ' num2str(opts.bb_crit) ]);
    disp(['     Broadband Passband: ' num2str(opts.bb_pass) ' Hz']);
    disp(['     Broadband Detrend: '  char(string(opts.hf_detrend))]);
    disp(' ')

    if opts.slope_test
        disp(['     Slope Test: slope > ' num2str(opts.slope_crit)]);
    else
        disp('     Slope Test: False');
    end
    disp('    ');
end

%% Get bad indicies
if opts.slope_test
    [spect, stimes, sfreqs] = multitaper_spectrogram_mex(data, Fs, [1 min(55,Fs/2)],[10,19],[10 5],[],[],[],false,false);
    B = zeros(length(stimes),2);
    for ii = 1:length(stimes)
        B(ii,:) = polyfit(log(sfreqs)',log(spect(:,ii)),1);
    end
    t = (0:length(data)-1)/Fs;
    bad_slope = interp1(stimes,double(B(:,1)>opts.slope_crit)',t,'nearest')';
    bad_slope(isnan(bad_slope)) = 1;
else
    bad_slope = zeros(size(data));
end

%Get bad indices
[~, ~, ~, is_flat] = get_chunks(data, Fs);
bad_inds = isnan(data) | isinf(data) | is_flat | bad_slope;
bad_inds = find_outlier_noise(data, bad_inds);

%Interpolate big gaps in data
t = 1:length(data);
good_data = data(~bad_inds);
data_fixed = interp1([0, t(~bad_inds), length(data)+1], [good_data(1); good_data; good_data(end)], t)';

%% Get high frequency artifacts
hf_artifacts = compute_artifacts(opts.hpFilt_high, opts.detrend_duration, opts.hf_crit, data_fixed, opts.smooth_duration, Fs, bad_inds, opts.isexcluded, opts.verbose,...
    'high frequency', opts.histogram_plot, opts.zscore_method, opts.hf_detrend, opts.diagnostic_plot);

%% Get broad band frequency artifacts
bb_artifacts = compute_artifacts(opts.hpFilt_broad, opts.detrend_duration, opts.bb_crit, data_fixed, opts.smooth_duration, Fs, bad_inds, opts.isexcluded, opts.verbose,...
    'broadband frequency', opts.histogram_plot, opts.zscore_method, opts.bb_detrend, opts.diagnostic_plot);

%% Join artifacts from different frequency bands and put back to original size
% artifacts = true(size(data));
artifacts = hf_artifacts | bb_artifacts | bad_inds;

%Block out artifacts if requested for exclusion
if strcmpi(opts.exclude_mode,'artifact')
    artifacts = artifacts | opts.isexcluded;
end

%% Add buffer on both sides of detected artifacts
if opts.buffer_duration > 0
    [cons, inds] = consecutive_runs(artifacts);
    for ii = 1:length(cons)
        buffer_start_idx = ceil(max(1, inds{ii}(1)-opts.buffer_duration*Fs));
        buffer_end_idx = ceil(min(length(artifacts), inds{ii}(end)+opts.buffer_duration*Fs));
        artifacts(buffer_start_idx:buffer_end_idx) = true;
    end
end

%% Sanity check before outputting
assert(length(artifacts) == length(data), 'Data vector length is inconsistent. Please check.')

%Force artifacts to be a row vector
if iscolumn(artifacts)
    artifacts = transpose(artifacts);
end

if opts.verbose
    disp(['Total time: ' num2str(toc(t_start)) ' seconds']);
    disp(' ');
end
end


function bad_inds = find_outlier_noise(data, bad_inds)
%% %Find all time points that have outlier high or low noise
data_mean = mean(data(~bad_inds));
data_std = std(data(~bad_inds));

outlier_scalar = 10;
low_thresh = data_mean - outlier_scalar * data_std;
high_thresh = data_mean + outlier_scalar * data_std;

inds = data <= low_thresh | data >= high_thresh;
bad_inds(inds) = true;
end


function detected_artifacts = compute_artifacts(filter_coeff, detrend_duration, crit, data_fixed,...
    smooth_duration, Fs, bad_inds, isexcluded, verbose, verbosestring, histogram_plot, zscore_method, detrend_on, diagnostic_plot)
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

%Perform a zero-phase filter
y_signal = filtfilt(filter_coeff, data_fixed);

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

    if diagnostic_plot
        axes(ax(5))
        plot(t,y_signal)
        title('Detrended Log Smoothed Magnitude');
    end
end

% create an indexing array marking bad timepoints before z-scoring
detected_artifacts = bad_inds;

%Take z-score
ysig = y_signal(~detected_artifacts & ~isexcluded);
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
over_crit = abs(y_signal)>crit & ~detected_artifacts & ~isexcluded;

%Loop until nothing over criterion
count = 1;
while any(over_crit)

    %Update the detected artifact time points
    detected_artifacts(over_crit & ~bad_inds) = true;

    %Compute modified z-score
    ysig = y_signal(~detected_artifacts & ~isexcluded);
    if strcmpi(zscore_method,'standard')
        ymid = mean(ysig);
        ystd = std(ysig);
    else
        ymid = median(ysig);
        ystd = mad(ysig);
    end

    y_signal = (y_signal - ymid)/ystd;

    %Find new criterion
    over_crit = abs(y_signal)>crit & ~detected_artifacts & ~isexcluded;

    if histogram_plot
        axes(ah); %#ok<LAXES>
        histogram(y_signal(~detected_artifacts & ~isexcluded), 100);
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
end
