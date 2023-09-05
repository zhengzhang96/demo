function [fh_baseband,Data_Matrix_1,bin_length,range_vector,total_time,info_1,frame_stamp_vec,queue,time_toc,recording_actual_time,frame_area_tot,range_vector_length]=rec_and_vis_data_1Radar_setDuration_IQ(xep_1,frame_start,frame_stop,actual_fps,fps,total_duration,visualise_data,visualise_fps,calibration,freqRes_STFT)

% Initial settings - Read and plot data from modules. 

% Settings for radar history plot
dt = 1/actual_fps; %Time for each frame in seconds.
history_duration=2; %In seconds
dt_plot = 50*dt;
plot_steps = floor(dt_plot*actual_fps);

if ~calibration
    total_duration = total_duration + history_duration;
end

% Compute steps in samples to visualise data (or update waitbar) at the
% fixed fps defined by visualise_fps.
if visualise_fps<actual_fps
    visualise_steps = ceil(actual_fps/visualise_fps); %in samples
else % in case the actual fps or recording is lower than the visualisation fps.
    error('Visualise frame-rate must be lower than recording frame-rate.')
end

% Initialise figure handle as empty 
fh_baseband = []; % In case we do not choose to visualise data, this will remain empty.

if ~visualise_data
    if ~calibration
        h = waitbar(0,'Recording Data...');
    else
        h = waitbar(0,'Recording Empty Room Data...');
    end
end

% Generate range vector
bin_length = 8 * 1.5e8/23.328e9; % range_decimation_factor * (c/2) / fs.
range_vector = (frame_start-1e-5):bin_length:(frame_stop+1e-5); % +-1e-5 to account for float precision.

frame_area_tot = (frame_stop+1e-5)-(frame_start-1e-5);
range_vector_length = length(range_vector);

%Data Matrix Initialisation
tot_frame_num = ceil(total_duration*actual_fps);
Data_Matrix_1 = single(NaN(tot_frame_num,length(range_vector)));

if ~calibration
    Data_Matrix_1_unfilt = single(NaN(tot_frame_num,length(range_vector)));
    Data_Matrix_1_plot = single(zeros(tot_frame_num,length(range_vector)));
    % Data_Matrix_1_plot = gpuArray(single(zeros(tot_frame_num,length(range_vector))));
end

% Initialise matrix to store frame stamps ( if frame is dropped, then
% info_1 will be zero anyway.
frame_stamp_vec = zeros(1,tot_frame_num);

% Creating a running-time variable
total_time = 0; %Setting timer to zero

if ~calibration
    % Settings for Doppler Spectrum plot
    if (rem(actual_fps,freqRes_STFT)/freqRes_STFT)>=0.5 || (rem(actual_fps,freqRes_STFT)/freqRes_STFT)==0
        f_num_STFT = floor(actual_fps/freqRes_STFT)+1;
    else
        f_num_STFT = floor(actual_fps/freqRes_STFT);
    end
    win_size = 100;
    window = kaiser(win_size,15);
    window_mtx = repmat(window,1,size(Data_Matrix_1,2)); % Arrange window as a matrix to be applied to all range bins simultaneously.
    noverlap = 75;
    % determine frame shift (step) in samples
    step = win_size - noverlap;
    % win_centre_ind = floor(win_size/2)+1 : step : size(Data_Matrix_1,1)-(floor(win_size/2))+1; % Corresponding to the center of the window
    win_upper_ind = win_size:step:size(Data_Matrix_1,1);
    win_lower_ind = win_upper_ind-win_size+1;

    f_STFT = (-actual_fps/2):actual_fps/f_num_STFT: (actual_fps/2) - (actual_fps/f_num_STFT);
    f_c = 7.29e9; %center frequency of carrier sinusoid.
    c=3e8; %Speed of light in m/s.
    lambda = c/f_c;
    v_STFT = (lambda/2)*f_STFT;
    % t_STFT = win_upper_ind*dt;
end

% Note: To find screen size use:
% set(0,'units','pixels')  
% %Obtains this pixel information
% Pix_SS = get(0,'screensize')

% If we selected to visualise data - Faster to pre-initialise figures
if visualise_data
    fh_baseband = figure('Name','Recording Baseband signals','NumberTitle','off','Position', [218 62 1100 700]);
    clf(fh_baseband); %clf(fig) clears the single figure with handle fig.

    % Create handles for update of plots.
    % Plotting the received signal
    subplot(4,4,[1,2,3,4])
    ph_signal = plot(NaN,NaN); %ph stands for plot handle
    th_signal = title(''); %th stands for title handle
    th_signal.String = 'Received Signal';
    ax_signal = gca; %ax stands for axes handle
    ylim(ax_signal,[0 0.03]);
    xlim([range_vector(1) range_vector(end)]);
    grid on;
    xlabel('Range [m]');
    ylabel('Magnitude (AU)')
    ph_signal.XData = range_vector; %x-axis values
    ph_signal.YData = Data_Matrix_1(1,:)';

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
    xlim([range_vector(1) range_vector(end)]);
    ylim([f_STFT(1) f_STFT(end)]);
    ph_doppler.XData = range_vector; %x-axis values
    ph_doppler.YData = f_STFT; %y-axis values
    ph_doppler.CData = 20*log10(eps*ones(length(f_STFT),length(range_vector))); %y-axis values
    
    yyaxis right % Showing y-axis in m/s
    ph_doppler_right = imagesc(0,0,0);
    ax_doppler_right = gca;
    set(ax_doppler_right,'YDir','normal')
    ph_doppler_right.YData = v_STFT; %y-axis values
    ph_doppler_right.XData = range_vector; %x-axis values
    xlim([range_vector(1) range_vector(end)]);
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
    th_history.String = 'Received Signal History';
    ax_history = gca;
    set(ax_history,'YDir','normal')
    set(ax_history,'CLim',[0 0.1]); % Set color limits of pcolor
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
    cal_filename = '../../Matlab Codes/Recording Codes/Recordings MAT files/Radar_Data_Demo_IQ_Calibration.mat'; % Same filename always
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

if ~calibration
    % Construct a 2nd order Butterworth High pass filter
    v_c = 0.2; %In m/s
    f_c = (2*v_c)/lambda; % Cut-off frequency in Hz.
    [b,a] = butter(2,f_c/(actual_fps/2),'high');
end

% Store number of frames in buffer - For checking recording quality.
queue=zeros(1,size(Data_Matrix_1,1));

% % Consider matched filtering
% % Transmitted pulse template - settings
% BW = 1.4e9; %Bandwidth (-10dB bandwidth more precisely)
% f_B = BW/2;
% tau = 1/(2*pi*f_B*sqrt(log10(exp(1))));
% V_tx = 4; %The value of this does not matter in the cross-correlation output.
% pulse_start_pos = 2*((0.18)+(bin_length/2))/c;
% 
% %Total length of pulse (i.e. PRI)
% PRF = 15.1875e6;
% R_max = (c/2)*(1/PRF);
% 
% % radar parameters
% dec_fact = 8; %Make sure it matches the one in RadarReturnsFromWalkingHuman_noR2_Scaling_BodyPart_Selection() function.
% fs_fast = 23.328e9/dec_fact; %Fast time sampling rate for X4 Radar
% 
% % Fast-time axis
% r_fast = (0:bin_length:(R_max-bin_length))-(bin_length/2);
% t_fast = (2/c)*r_fast;
% [~,ind_start] = min(abs(r_fast-frame_start));
% [~,ind_end] = min(abs(r_fast-frame_stop));
% 
% % Construct fast-time frequency axis - Using zero-padding
% extra_precision = 1; %integer value.
% nFFT_fast = 2^(nextpow2(length(r_fast))+extra_precision); %Total length with zero-padding.
% %frequency axis
% fscale_fast=((-fs_fast/2):(fs_fast/nFFT_fast):(fs_fast/2-(fs_fast/nFFT_fast))).';
% 
% %Obtaining frequency-shifted Gaussian pulse - Our template
% gaussian_pulse = V_tx * exp(-((t_fast-pulse_start_pos).^2)./(2*(tau^2)));

% Store time needed per loop of code.
time_toc = zeros(1,tot_frame_num);

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

while (total_time<total_duration)

    tic
    [~,data_length_1,info_1,data_1] = xep_1.read_message_data_float();
    info_1 = info_1 - info_1_old;
    queue(info_1)= xep_1.peek_message_data_float();

    % Generate IQ vector - for baseband data.
    i_vec_1 = data_1(1:data_length_1/2);
    q_vec_1 = data_1(data_length_1/2+1:data_length_1);
    iq_vec_1 = i_vec_1 + 1i*q_vec_1;

    % Remove empty room average vector when not performing empty room
    % recordings.  Also remove online average vector (DC).
    if ~calibration
        iq_vec_1 = iq_vec_1 - (calibration_vector.'); %Empty room subtraction
    end

    %     % Apply matched filtering based on transmitted pulse template defined
    %     % above
    %     % Add zeros to range bins that were not included in the recorded data (to
    %     % match the size of the pulse template)
    %     single_frame_temp = [zeros(ind_start-1,1);iq_vec_1;zeros(length(r_fast)-ind_end,1)];
    %
    %     [rxy,lags] = xcorr(single_frame_temp,gaussian_pulse.');
    %     matched_filter_applied = rxy(lags>=0);
    %
    %     % Dealing with the time shift corresponding to 0.18m
    %     matched_filter_applied_fft = fftshift(fft(matched_filter_applied,nFFT_fast)).*exp(-1i*2*pi*fscale_fast*pulse_start_pos);
    %     matched_filter_applied_temp = ifft(ifftshift(matched_filter_applied_fft));
    %     matched_filter_applied_shift = matched_filter_applied_temp(1:length(matched_filter_applied));
    %
    %     iq_vec_1 = matched_filter_applied_shift(ind_start:ind_end);

    if ~calibration
        % Store each radar frame in Data Matrix - we store the iq_vectors
        % rows = fast-time (range) dimension
        % columns = slow-time (frame) dimension
        Data_Matrix_1_unfilt(info_1,:) = iq_vec_1.'; %Attach each vector (iq_vec' is row vector) in next row.

        % High Pass filtering stage - 2nd Order Butterworth
        if info_1>=3 % Order+1
            Data_Matrix_1(info_1,:) = b(1)*Data_Matrix_1_unfilt(info_1,:) + b(2)*Data_Matrix_1_unfilt(info_1-1,:) + b(3)*Data_Matrix_1_unfilt(info_1-2,:)...
                -a(2)*Data_Matrix_1(info_1-1,:) - a(3)*Data_Matrix_1(info_1-2,:);
        else
            Data_Matrix_1(info_1,:) = Data_Matrix_1_unfilt(info_1,:);
        end

        % Matrix for history data - Ensure that before history_duration seconds
        % we do not have any values since we do not plot these.
        if info_1>history_duration/dt
            Data_Matrix_1_plot(info_1,:)=abs(Data_Matrix_1(info_1,:));
        end
    else
        % Store each radar frame in Data Matrix - we store the iq_vectors
        % rows = fast-time (range) dimension
        % columns = slow-time (frame) dimension
        Data_Matrix_1(info_1,:)=iq_vec_1.'; %Attach each vector (iq_vec' is row vector) in next row.
    end


    %     % 3-pulse canceller
    %     if info_1>=3
    %         Data_Matrix_1(info_1,:) = Data_Matrix_1_unfilt(info_1,:) - 2* Data_Matrix_1_unfilt(info_1-1,:) + Data_Matrix_1_unfilt(info_1-2,:);
    %     end
    %     Data_Matrix_1_plot(info_1,:)=abs(Data_Matrix_1(info_1,:));

    %     % 2-pulse canceller
    %     if info_1>=2
    %         Data_Matrix_1(info_1,:) = Data_Matrix_1_unfilt(info_1,:) - Data_Matrix_1_unfilt(info_1-1,:);
    %     end
    %     Data_Matrix_1_plot(info_1,:)=abs(Data_Matrix_1(info_1,:));

    if ~calibration
        % Apply STFT in real-time - FFT on data*window function.
        [TF,loc]=ismember(info_1,win_upper_ind);
        if TF && ~mod(info_1,visualise_steps)
            % Multiply segment of signal by window function
            Data_Matrix_1_window = Data_Matrix_1(win_lower_ind(loc):win_upper_ind(loc),:).*window_mtx;
            % Compute the Fourier Transform of this segment of data
            STFT = fftshift(fft(Data_Matrix_1_window,f_num_STFT,1),1); % fft along first dimension.
        end
    end
    
    % Store frame stamps
    frame_stamp_vec(1,info_1)=info_1;

    %Update timer
    total_time=dt*info_1;

    % Update title of plot during intiialisation
    if info_1<=round(history_duration/dt,0) && visualise_data
        perc_completed = (info_1/round(history_duration/dt,0))*100;
        
%         general_title.String=['Info_1 = ' num2str(info_1)];
        if perc_completed>=(counter_initialisation*20)
            general_title.String=['Initialising Recording (' num2str(counter_initialisation*20) '%) ...'];
            drawnow;
            counter_initialisation = counter_initialisation+1;
        end
    end


    if visualise_data
        if ~mod(info_1,visualise_steps) && info_1>=win_upper_ind(1) && info_1<=win_upper_ind(end) % Plot at the defined visualisation frame rate

            if info_1>round(history_duration/dt,0)
                % Update data in handles for plotting

                % Plot received signal
                ph_signal.YData = abs(Data_Matrix_1(info_1,:))';
                peak_val = max([0.03 max(abs(Data_Matrix_1(info_1,:))) + 0.01]);
                ylim(ax_signal,[0 peak_val]);

                % Plot doppler information
                ph_doppler.CData = 20*log10(abs(conj(STFT))+eps);

                % Plot received signal history
                if total_time<=history_duration
                    ph_history.CData = Data_Matrix_1_plot(1:plot_steps:info_1,:);
                elseif total_time>history_duration
                    ph_history.CData = Data_Matrix_1_plot((info_1-floor(history_duration/dt)):plot_steps:info_1-plot_steps,:);
                end
                %                 set(ax_history,'CLim',[0 peak_val+0.01]); % Set color limits of pcolor

                general_title.String = sprintf('Recording Data - Elapsed Time = %.2f s',round(total_time - history_duration,2));

                drawnow;
            end
        end
    else
        waitbar(total_time/total_duration,h)
    end

    time_toc(info_1)=toc;
end

recording_actual_time{2} = datestr(now,'HH:MM:SS:FFF');

if ~visualise_data
    waitbar(1,h)
    pause(0.1); % For stability
    close(h)
end

% Truncate data and information showing data after history_duration has
% passed. Before this, we do not visualise any data.
if ~calibration
    total_time = total_time-history_duration;
    info_1 = info_1 - round(history_duration/dt);
    Data_Matrix_1 = Data_Matrix_1(round(history_duration/dt)+1:end,:);
    frame_stamp_vec = frame_stamp_vec(round(history_duration/dt)+1:end) - round(history_duration/dt);
    % Make negative values equal to zero.
    frame_stamp_vec = max(frame_stamp_vec,0);
    queue = queue(round(history_duration/dt)+1:end);
    time_toc = time_toc(round(history_duration/dt)+1:end);

    % Shifting recording_actual_time{1} by history_duration seconds.
    old_start_time = datetime(recording_actual_time{1}(1:8),'InputFormat','HH:mm:ss'); %Convert to datetime format
    new_start_time = old_start_time + seconds(2);
    recording_actual_time{1} = append(datestr(new_start_time,'HH:MM:SS'),recording_actual_time{1}(9:end));
end


fprintf('\nTotal elapsed time = %f s\n\n',total_time);

end