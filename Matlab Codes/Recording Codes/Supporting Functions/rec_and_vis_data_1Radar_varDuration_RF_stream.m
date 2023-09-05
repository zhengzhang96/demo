function [fh_RF]=rec_and_vis_data_1Radar_varDuration_RF_stream(xep_1,frame_start,frame_stop,actual_fps,fps,visualise_data,visualise_fps,calibration,freqRes_STFT)

% Initial settings - Read and plot data from modules.

% Settings for radar history plot
dt = 1/actual_fps; %Time for each frame in seconds.
history_duration=2; %In seconds
dt_plot = 50*dt;
plot_steps = floor(dt_plot*actual_fps);

% Compute steps in samples to visualise data (or update waitbar) at the
% fixed fps defined by visualise_fps.
if visualise_fps<actual_fps
    visualise_steps = ceil(actual_fps/visualise_fps); %in samples
else % in case the actual fps or recording is lower than the visualisation fps.
    error('Visualise frame-rate must be lower than recording frame-rate.')
end

% Initialise figure handle as empty 
fh_RF = []; % In case we do not choose to visualise data, this will remain empty.

if ~visualise_data
    if ~calibration
        h = waitbar(0,'Recording Data...');
    else
        h = waitbar(0,'Recording Empty Room Data...');
    end
end

% Generate range vector - We need no decimation in RF signal.
bin_length = 1.5e8/23.328e9; % range_decimation_factor * (c/2) / fs.
range_vector = (frame_start-1e-5):bin_length:(frame_stop+1e-5)-bin_length; % +-1e-5 to account for float precision.

% Data Matrix Initialisation to store frames for HP filter
Data_Matrix_1_unfilt = single(NaN(3,length(range_vector))); % HP filter order is 3.
Data_Matrix_1 = single(NaN(3,length(range_vector)));

tot_frame_num_history = round(history_duration/dt);
Data_Matrix_1_plot = single(zeros(tot_frame_num_history,length(range_vector)));

% Initialise variable to store sum of dropped frames ( if frame is dropped, then
% info_1 will be zero).
frames_stamp_vec = zeros(1,round(60/dt)); % Give a check for dropped frames every 1 minute of recording.

% Store number of frames in buffer - For checking recording quality (every minute).
queue=zeros(1,round(60/dt));

% Creating a running-time variable
total_time = 0; %Setting timer to zero

f_0 = 7.29e9; %Center frequency
c=3e8; %Speed of light in m/s.

% Settings for RF to IQ conversion
t_fast_conversion = (2/c)*(range_vector)';
complex_fact = 2*exp(-2j*pi*f_0*t_fast_conversion);
dec_fact_conv = 4;
% Design the FIR LPF used in downconversion
% The FIR LPF is designed based on the method resample in MATLAB. 
% The ideal antialiasing filter has normalized cutoff frequency fc = 1/downsampling_factor (π rad/sample)
% and gain (p) of 1. 
p=1;
q=dec_fact_conv;
f_c_lp = 1/q;
order = 25;
b_lp = fir1(order,f_c_lp);
b_lp = p*b_lp/sum(b_lp);
a_lp=1;
%Initialise the data matrix to store the IQ data
range_vector_IQ = range_vector(1:q:end);
% Data_Matrix_1_IQ = single(NaN(tot_frame_num,length(range_vector_IQ)));

% Settings for Doppler Spectrum plot
if (rem(actual_fps,freqRes_STFT)/freqRes_STFT)>=0.5 || (rem(actual_fps,freqRes_STFT)/freqRes_STFT)==0
    f_num_STFT = floor(actual_fps/freqRes_STFT)+1;
else
    f_num_STFT = floor(actual_fps/freqRes_STFT);
end
win_size = 100;
window = kaiser(win_size,15);
window_mtx = repmat(window,1,length(range_vector_IQ)); % Arrange window as a matrix to be applied to all range bins simultaneously.
noverlap = 75;
% determine frame shift (step) in samples
step = win_size - noverlap;
win_upper_ind_start = win_size+tot_frame_num_history;
win_lower_ind_start = win_upper_ind_start-win_size+1;

upper_ind = win_upper_ind_start;

Data_Matrix_1_STFT = single(zeros(win_size,length(range_vector_IQ)));

f_STFT = (-actual_fps/2):actual_fps/f_num_STFT: (actual_fps/2) - (actual_fps/f_num_STFT);
lambda = c/f_0;
v_STFT = (lambda/2)*f_STFT;

% % Store STFT tensor for averaging every two frames along time.
STFT_tensor = zeros(2,f_num_STFT,length(range_vector_IQ));
STFT_counter = 1; % Counter for computer average of STFT tensor.

% Note: To find screen size use:
% set(0,'units','pixels')  
% %Obtains this pixel information
% Pix_SS = get(0,'screensize')

% If we selected to visualise data
if visualise_data
    fh_RF = figure('Name','Recording RF signals','NumberTitle','off','Position',[218 62 1100 700],'KeyPressFcn',@key_pressed_fcn_ENTER);
    clf(fh_RF); %clf(fig) clears the single figure with handle fig.

    % Create handles for update of plots.
    % Plotting the received signal
    subplot(4,4,[1,2,3,4])
    ph_signal_RF = plot(NaN,NaN); %ph stands for plot handle
    th_signal = title(''); %th stands for title handle
    th_signal.String = 'Received Signal ({\color{blue}RF},{\color{red}IQ})';
    ax_signal = gca; %ax stands for axes handle
    ylim(ax_signal,[-0.02 0.02]);
    xlim([0.2 range_vector(end)]);
    grid on;
    xlabel('Range [m]');
    ylabel('Amplitude (AU)')
    ph_signal_RF.XData = range_vector; %x-axis values
    ph_signal_RF.YData = zeros(1,length(range_vector))';
    hold on
    ph_signal_IQ = plot(NaN,NaN); %ph stands for plot handle
    ph_signal_IQ.XData = range_vector_IQ; %x-axis values
    ph_signal_IQ.YData = zeros(1,length(range_vector_IQ))';
%     legend('RF','IQ') -> This causes code to slow down

    % Plotting the doppler information
    subplot(4,4,[7,8,11,12,15,16])
    yyaxis left % Showing y-axis in Hz
    ph_doppler = imagesc(0,0,0);
    th_doppler = title('');
    th_doppler.String = 'Range-Doppler Information';
    ax_doppler = gca;
    set(ax_doppler,'YDir','normal')
    set(ax_doppler,'CLim',[-35 5]);
    colormap(jet)
    xlabel('Range [m]'); ylabel('Doppler Frequency (Hz)'); 
    xlim([range_vector_IQ(1) range_vector_IQ(end)]);
    ylim([f_STFT(1) f_STFT(end)]);
    ph_doppler.XData = range_vector_IQ; %x-axis values
    ph_doppler.YData = f_STFT; %y-axis values
    ph_doppler.CData = 20*log10(eps*ones(length(f_STFT),length(range_vector_IQ))); %y-axis values
    
    yyaxis right % Showing y-axis in m/s
    ph_doppler_right = imagesc(0,0,0);
    ax_doppler_right = gca;
    set(ax_doppler_right,'YDir','normal')
    ph_doppler_right.YData = v_STFT; %y-axis values
    ph_doppler_right.XData = range_vector_IQ; %x-axis values
    xlim([range_vector_IQ(1) range_vector_IQ(end)]);
    ylim([v_STFT(1) v_STFT(end)]);
    ylabel('Velocity (m/s)');
    cla

    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';

    % Plotting the received signal history
    subplot(4,4,[5,6,9,10,13,14])
    ph_history = imagesc(0,0,0);
    th_history = title('');
    th_history.String = 'Received Signal History (Rectified)';
    ax_history = gca;
    set(ax_history,'YDir','normal')
    set(ax_history,'CLim',[0 0.08]); % Set color limits of pcolor
    colormap(jet)
    xlabel('Range (m)'); ylabel('Past Time (s)');
    ph_history.XData = range_vector; %x-axis values
    xlim([range_vector(1) range_vector(end)]);
    ph_history.YData = (-history_duration+dt_plot:dt_plot:0);
    ylim([-history_duration+dt_plot 0]);
    ph_history.CData = Data_Matrix_1_plot(1:plot_steps:round(history_duration/dt),:);

    general_title = sgtitle ('Initialising Recording (0%) ...');
    general_title.FontSize = 12;

end

% Consider empty room subtraction
if ~calibration %Only when not in calibration process.
    cal_filename = '../../Matlab Codes/Recording Codes/Empty Room Recordings MAT files/Radar_Data_RF_Calibration.mat'; % Same filename always
    if isfile(cal_filename)
        calibration_matrix = open(cal_filename).Data_Matrix_1;
        calibration_range_vector = open(cal_filename).range_axis;
        if isequal(calibration_range_vector,single(range_vector)) % Data was stored as 'single'
            calibration_vector = mean(calibration_matrix,1,'omitnan');
        else
            calibration_vector = zeros(1,length(range_vector));
        end
    else
        calibration_vector = zeros(1,length(range_vector));
    end
end

% Construct a 2nd order Butterworth High pass filter
v_c_hp = 0.2; %In m/s
f_c_hp = (2*v_c_hp)/lambda; % Cut-off frequency in Hz.
[b_hp,a_hp] = butter(2,f_c_hp/(actual_fps/2),'high');

% Construct template for matched filtering
% Transmitted pulse template - settings
BW = 1.4e9; %Bandwidth (-10dB bandwidth more precisely)
f_B = BW/2;
tau = 1/(2*pi*f_B*sqrt(log10(exp(1))));
V_tx = 0.08; %The value of this does not matter in the cross-correlation output.
pulse_start_pos = 2*((0.18))/c;

%Total length of pulse (i.e. PRI)
PRF = 15.1875e6;
R_max = (c/2)*(1/PRF);

% radar parameters
dec_fact_matched_filt = 1; %Make sure it matches the one in RadarReturnsFromWalkingHuman_noR2_Scaling_BodyPart_Selection() function.
fs_fast = 23.328e9/dec_fact_matched_filt; %Fast time sampling rate for X4 Radar

% Fast-time axis
r_fast = 0:bin_length:(R_max-bin_length);
t_fast = (2/c)*r_fast;
[~,ind_start] = min(abs(r_fast-frame_start));
[~,ind_end] = min(abs(r_fast-frame_stop));

% Construct fast-time frequency axis - Using zero-padding
extra_precision = 1; %integer value.
nFFT_fast = 2^(nextpow2(length(r_fast))+extra_precision); %Total length with zero-padding.
%frequency axis
fscale_fast=((-fs_fast/2):(fs_fast/nFFT_fast):(fs_fast/2-(fs_fast/nFFT_fast))).';

%Obtaining frequency-shifted Gaussian pulse - Our template
gaussian_pulse = V_tx * exp(-((t_fast-pulse_start_pos).^2)./(2*(tau^2)));
gaussian_pulse_shifted = gaussian_pulse.*cos(2*pi*f_0*(t_fast-pulse_start_pos));

if visualise_data
    drawnow;
end

xep_1.x4driver_set_fps(0);
% xep_1.x4driver_get_fps() % In case you need to check
% Clear message buffers
while xep_1.peek_message_data_float > 0
    [~,~,info_1_old,~] = xep_1.read_message_data_float();
end
xep_1.x4driver_set_fps(fps);
% xep_1.x4driver_get_fps() % In case you need to check

% Loop until Total Duration set by the User has elapsed.
recording_actual_time{1} = datestr(now,'HH:MM:SS:FFF');

% Counter for assisting in printing initialisation completion.
counter_initialisation = 1;

global keep %#ok<GVMIS> 
keep = 1; % Variable will become 0 if we press any key.

disp([newline 'After initialization, press ENTER (RETURN) key to stop recording...' newline])

while keep
    
    [~,~,info_1,data_1] = xep_1.read_message_data_float();
    info_1 = info_1 - info_1_old;

    % Remove empty room average vector when not performing empty room
    % recordings.  Also remove online average vector (DC).
    if ~calibration
        data_1 = data_1 - (calibration_vector.'); %Empty room subtraction
    end
    
    % Apply matched filtering based on transmitted pulse template defined
    % above
    % Add zeros to range bins that were not included in the recorded data (to
    % match the size of the pulse template)
    single_frame_temp = [zeros(ind_start,1);data_1;zeros(length(r_fast)-ind_end,1)];
    [rxy,lags] = xcorr(single_frame_temp,gaussian_pulse_shifted.');
    matched_filter_applied = rxy(lags>=0);

    % Dealing with the time shift corresponding to 0.18m
    matched_filter_applied_fft = fftshift(fft(matched_filter_applied,nFFT_fast)).*exp(-1i*2*pi*fscale_fast*pulse_start_pos);
    matched_filter_applied_temp = ifft(ifftshift(matched_filter_applied_fft));
    matched_filter_applied_shift = matched_filter_applied_temp(1:length(matched_filter_applied));
    
    data_1 = real(matched_filter_applied_shift(ind_start+1:ind_end));

    %Store each radar frame in Data Matrix
    % rows = slow-time (frame) dimension
    % columns = fast-time (range) dimension
    % Apply circular shift to the stored frames, so that newest frame is at
    % position with index 3. (Current and two past frames)
    Data_Matrix_1_unfilt = circshift(Data_Matrix_1_unfilt,-1,1);
    Data_Matrix_1_unfilt(3,:) = data_1'; %Attach each vector (data' is row vector) in next row.

    % High Pass filtering stage - 2nd Order Butterworth
    Data_Matrix_1 = circshift(Data_Matrix_1,-1,1);
    if info_1>=3 % Order+1
        Data_Matrix_1(3,:) = b_hp(1)*Data_Matrix_1_unfilt(3,:) + b_hp(2)*Data_Matrix_1_unfilt(2,:) + b_hp(3)*Data_Matrix_1_unfilt(1,:)...
            -a_hp(2)*Data_Matrix_1(2,:) - a_hp(3)*Data_Matrix_1(1,:);
    else
        Data_Matrix_1(3,:) = Data_Matrix_1_unfilt(3,:);
    end
    
    % Matrix for history data - Ensure that before history_duration seconds
    % we do not have any values since we do not plot these.
    if info_1>history_duration/dt
        Data_Matrix_1_plot = circshift(Data_Matrix_1_plot,-1,1);
        Data_Matrix_1_plot(end,:)=Data_Matrix_1(3,:);
    end

    if info_1>=win_lower_ind_start
        % Convert data from IQ to RF
        RF_data_complex = (Data_Matrix_1(3,:))'.*complex_fact;
        % Apply filter - Note: filtfilt operates along the first array dimension with size greater than 1.
        RF_data_complex_filt = FiltFiltM(b_lp,a_lp,RF_data_complex);
        IQ_data = conj(RF_data_complex_filt(1:q:end,:)');
        Data_Matrix_1_STFT = circshift(Data_Matrix_1_STFT,-1,1);
        Data_Matrix_1_STFT(end,:) = IQ_data;
    end

    % Apply STFT in real-time - FFT on data*window function.
    if (info_1==upper_ind) %&& ~mod(info_1,visualise_steps)
        % Multiply segment of signal by window function
        Data_Matrix_1_window = Data_Matrix_1_STFT.*window_mtx;
        % Compute the Fourier Transform of this segment of data
        STFT = fftshift(fft(Data_Matrix_1_window,f_num_STFT,1),1); % fft along first dimension.
        STFT_tensor = circshift(STFT_tensor,-1,1);
        STFT_tensor(2,:,:) = STFT;
        if STFT_counter>=2
            STFT = squeeze(mean(STFT_tensor,1));
        end
        STFT_counter= STFT_counter+1;
        upper_ind = upper_ind + step;
    end

    % Store data stamps
    frames_stamp_vec = circshift(frames_stamp_vec,-1);
    frames_stamp_vec(end) = info_1;

    % Store number of frames in buffer
    queue = circshift(queue,-1);
    queue(end) = xep_1.peek_message_data_float();

    % Check how many frames are dropped every 1 minute
    if ~mod(info_1,round(60/dt))
        frames_dropped_sum = sum(frames_stamp_vec==0);
        disp(['# Frames dropped from Radar in last minute: ' num2str(frames_dropped_sum) ' (' num2str((frames_dropped_sum/round(60/dt))*100) '%)']);
        disp(['Mean/Std of Frames in Buffer in last minute: ' num2str(mean(queue)) ' (' num2str(mean(queue)*dt*1000) ' ms) / ' num2str(std(queue)) ' (' num2str(std(queue)*dt*1000) ' ms)']);
        disp(['Max/Min of Frames in Buffer in last minute: ' num2str(max(queue)) ' (' num2str(max(queue)*dt*1000) ' ms) / ' num2str(min(queue)) ' (' num2str(min(queue)*dt*1000) ' ms)' newline]);
    end
    
    %Update timer
    total_time=dt*info_1;

    % Update title of plot during intiialisation
    if info_1<=round(history_duration/dt,0)
        perc_completed = (info_1/round(history_duration/dt,0))*100;

        %         general_title.String=['Info_1 = ' num2str(info_1)];
        if perc_completed>=(counter_initialisation*20)
            general_title.String=['Initialising Recording (' num2str(counter_initialisation*20) '%) ...'];
            drawnow;
            counter_initialisation = counter_initialisation+1;
        end
    end
    
    if ~mod(info_1,visualise_steps) && info_1>=win_upper_ind_start % Plot at the defined visualisation frame rate
        if visualise_data

            if info_1>round(history_duration/dt,0)
                % Update data in handles for plotting
                
                % Plot received signal 
                ph_signal_RF.YData = Data_Matrix_1(3,:)';
                ph_signal_IQ.YData = abs(IQ_data)';
                peak_val = max([0.03 max(abs(Data_Matrix_1(3,:))) + 0.01]);
                ylim(ax_signal,[min([-0.02 (-peak_val - 0.01)]) max([0.02 (peak_val + 0.01)])]);

                % Plot doppler information
                ph_doppler.CData = 20*log10(abs(conj(STFT))+eps);

                % Plot received signal history
                ph_history.CData = abs(Data_Matrix_1_plot(1:plot_steps:end,:));

                general_title.String = sprintf(['{\\bfRecording Radar Data}' newline 'Date: ' datestr(now, 'dd/mm/yyyy'),' , Time: ',datestr(now,'HH:MM:SS:FFF')]);

                drawnow;
            end
        else
            waitbar(total_time/total_duration,h)
        end
    end
end

recording_actual_time{2} = datestr(now,'HH:MM:SS:FFF'); %#ok<*UNRCH> 

if ~visualise_data
    waitbar(1,h)
    pause(0.01); % For stability
    close(h)
end

total_time = total_time-history_duration;
old_start_time = datetime(recording_actual_time{1}(1:8),'InputFormat','HH:mm:ss'); %Convert to datetime format
new_start_time = old_start_time + seconds(2);
recording_actual_time{1} = append(datestr(new_start_time,'HH:MM:SS'),recording_actual_time{1}(9:end));

fprintf('\nTotal elapsed time = %.1f s\n\n',total_time);

disp(['Actual Recording Start Time:  ' recording_actual_time{1}]);
disp(['Actual Recording Stop Time:  ' recording_actual_time{2} newline]);

end