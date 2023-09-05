function [fh_baseband]=rec_and_vis_data_1Radar_varDuration_IQ_stream_working_3(xep_1,frame_start,frame_stop,actual_fps,fps,visualise_data,visualise_fps,calibration,freqRes_STFT)

% Initial settings - Read and plot data from module.

% Settings for radar history plot
dt = 1/actual_fps; %Time for each frame in seconds.
history_duration=10; %In seconds
dt_plot = 50*dt;
plot_steps = floor(dt_plot*actual_fps);

% Compute steps in samples to visualise data (or update waitbar) at the
% fixed fps defined by visualise_fps.
if visualise_fps<actual_fps
    visualise_steps = ceil(actual_fps/visualise_fps); %in samples
else
    error('Visualise frame-rate must be lower than recording frame-rate.')
end

% Initialise figure handle as empty 
fh_baseband = []; % In case we do not choose to visualise data, this will remain empty.

if ~visualise_data
    h = waitbar(0,'Recording Empty Room Data...');
end

% Generate range vector
bin_length = 8 * 1.5e8/23.328e9; % range_decimation_factor * (c/2) / fs.
range_vector = (frame_start-1e-5):bin_length:(frame_stop+1e-5); % +-1e-5 to account for float precision.

range_vector_length = length(range_vector);

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

% Settings for Doppler Spectrum plot
if (rem(actual_fps,freqRes_STFT)/freqRes_STFT)>=0.5 || (rem(actual_fps,freqRes_STFT)/freqRes_STFT)==0
    f_num_STFT = floor(actual_fps/freqRes_STFT)+1;
else
    f_num_STFT = floor(actual_fps/freqRes_STFT);
end
win_size = 100;
window = kaiser(win_size,15);
window_mtx = repmat(window,1,range_vector_length); % Arrange window as a matrix to be applied to all range bins simultaneously.
noverlap = 75;
% determine frame shift (step) in samples
step = win_size - noverlap;
win_upper_ind_start = win_size+tot_frame_num_history;
win_lower_ind_start = win_upper_ind_start-win_size+1;

upper_ind = win_upper_ind_start;

Data_Matrix_1_STFT = single(zeros(win_size,length(range_vector)));

f_STFT = (-actual_fps/2):actual_fps/f_num_STFT: (actual_fps/2) - (actual_fps/f_num_STFT);
f_c = 7.29e9; %center frequency of carrier sinusoid.
c=3e8; %Speed of light in m/s.
lambda = c/f_c;
v_STFT = (lambda/2)*f_STFT;

% % Store STFT tensor for averaging every two frames along time.
STFT_tensor = zeros(2,f_num_STFT,range_vector_length);
STFT_counter = 1; % Counter for computer average of STFT tensor.

r=1; %Contrast enhancement parameter

% Note: To find screen size use:
% set(0,'units','pixels')  
% %Obtains this pixel information
% Pix_SS = get(0,'screensize')

% If we selected to visualise data - Faster to pre-initialise figures
if visualise_data
    fh_baseband = figure('Name','Recording Baseband signals','NumberTitle','off','Position', [218 62 1100 700],'KeyPressFcn',@key_pressed_fcn_ENTER);
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
    ph_signal.YData = zeros(1,range_vector_length)';

    % Plotting the doppler information
    subplot(4,4,[7,8,11,12,15,16])
    yyaxis left % Showing y-axis in Hz
    ph_doppler = imagesc(0,0,0);
    th_doppler = title('');
    th_doppler.String = 'Range-Doppler Information';
    ax_doppler = gca;
    set(ax_doppler,'YDir','normal')
%     set(ax_doppler,'CLim',[-35 5]);
    set(ax_doppler,'CLim',[0.85 1]);
    colormap(jet)
    xlabel('Range [m]'); ylabel('Doppler Frequency (Hz)'); 
    xlim([range_vector(1) range_vector(end)]);
    ylim([f_STFT(1) f_STFT(end)]);
    ph_doppler.XData = range_vector; %x-axis values
    ph_doppler.YData = f_STFT; %y-axis values
    ph_doppler.CData = 20*log10(eps*ones(length(f_STFT),length(range_vector))); %y-axis values
    hold on
    ph_upper_envelope = plot(NaN,NaN,'k-','Linewidth',2.5);
    ph_upper_envelope.XData = range_vector;
    ph_upper_envelope.YData = zeros(1,length(range_vector));
    ph_lower_envelope = plot(NaN,NaN,'k-','Linewidth',2.5);
    ph_lower_envelope.XData = range_vector;
    ph_lower_envelope.YData = zeros(1,length(range_vector));
    ph_torso_loc = plot(NaN,NaN,'.','MarkerSize',20);
    
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
    set(ax_history,'CLim',[0 0.01]); % Set color limits of pcolor
    colormap(jet)
    xlabel('Range (m)'); ylabel('Past Time (s)');
    ph_history.XData = range_vector; %x-axis values
    xlim([range_vector(1) range_vector(end)]);
    ph_history.YData = (-history_duration+dt_plot:dt_plot:0);
    ylim([-history_duration+dt_plot 0]);
    ph_history.CData = Data_Matrix_1_plot(1:plot_steps:round(history_duration/dt),:);

    general_title = sgtitle ('Initialising Recording (0%) ...');
    general_title.FontSize = 14;
end

% Consider empty room subtraction
if ~calibration %Only when not in calibration process.
    cal_filename = '../../Matlab Codes/Recording Codes/Empty Room Recordings MAT files/Radar_Data_IQ_Calibration.mat'; % Same filename always
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
v_c = 0.2; %In m/s
f_c = (2*v_c)/lambda; % Cut-off frequency in Hz.
[b,a] = butter(2,f_c/(actual_fps/2),'high');

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

% Loop until Total Duration set by the User has elapsed.
recording_actual_time{1} = datestr(now,'HH:MM:SS:FFF');

% Counter for assisting in printing initialisation completion.
counter_initialisation = 1;

walking_speed = 0;
counter_walking = 0;

global keep %#ok<GVMIS> 
keep = 1; % Variable will become 0 if we press the ENTER (return) key.

disp([newline 'After initialization, press ENTER (RETURN) key to stop recording...' newline])

reset_steps = 1;

while keep
    
    [~,data_length_1,info_1,data_1] = xep_1.read_message_data_float();
    info_1 = info_1 - info_1_old;
    
    % Generate IQ vector - for baseband data.
    i_vec_1 = data_1(1:data_length_1/2);
    q_vec_1 = data_1(data_length_1/2+1:data_length_1);
    iq_vec_1 = i_vec_1 + 1i*q_vec_1;

    % Remove empty room average vector when not performing empty room
    % recordings.  Also remove online average vector (DC).
    if ~calibration
        iq_vec_1 = iq_vec_1 - (calibration_vector.'); %Empty room subtraction
    end
    
    % Store each radar frame in Data Matrix - we store the iq_vectors
    % rows = fast-time (range) dimension
    % columns = slow-time (frame) dimension
    % Apply circular shift to the stored frames, so that newest frame is at
    % position with index 3. (Current and two past frames)
    Data_Matrix_1_unfilt = circshift(Data_Matrix_1_unfilt,-1,1);
    Data_Matrix_1_unfilt(3,:) = iq_vec_1.'; %Attach each vector (iq_vec' is row vector) in next row.

    % High Pass filtering stage - 2nd Order Butterworth
    Data_Matrix_1 = circshift(Data_Matrix_1,-1,1);
    if info_1>=3 % Order+1 to avaoid calculations involving NaN's
        Data_Matrix_1(3,:) = b(1)*Data_Matrix_1_unfilt(3,:) + b(2)*Data_Matrix_1_unfilt(2,:) + b(3)*Data_Matrix_1_unfilt(1,:)...
            -a(2)*Data_Matrix_1(2,:) - a(3)*Data_Matrix_1(1,:);
    else
        Data_Matrix_1(3,:) = Data_Matrix_1_unfilt(3,:);
    end

    % Matrix for history data - Ensure that before history_duration seconds
    % we do not have any values since we do not plot these.
    if info_1>history_duration/dt
        Data_Matrix_1_plot = circshift(Data_Matrix_1_plot,-1,1);
        Data_Matrix_1_plot(end,:)=abs(Data_Matrix_1(3,:));
    end
    
    if info_1>=win_lower_ind_start
        Data_Matrix_1_STFT = circshift(Data_Matrix_1_STFT,-1,1);
        Data_Matrix_1_STFT(end,:) = Data_Matrix_1(3,:);
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
        %Contrast enhancement step
        mu = mean(abs(conj(STFT)),'all');
        std_val_new = std(abs(conj(STFT)),[],'all','omitnan');
        test = mu+3*std_val_new;
        if test<1e-3
            mu=1e-3;
        end
        STFT_new = (abs(conj(STFT)).^r)./((abs(conj(STFT)).^r) + (mu.^r));
        STFT_new_thresholded = STFT_new;
        STFT_new_thresholded(STFT_new_thresholded <= 0.85)=0;
        

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

                % Find torso in Doppler
                [max_all_range,range_indices] = max(STFT_new,[],2);
                [Torso_values,f_indices] = maxk(max_all_range,5);
                Torso_val = mean(Torso_values);
                f_ind = round(mean(f_indices));
                range_loc_torso = mean(range_vector(range_indices(f_indices)));

                ph_torso_loc.XData=NaN; ph_torso_loc.YData=NaN;

                % Extract upper and lower envelopes of each STFT frame
                upper_envelope = zeros(1,length(range_vector));
                lower_envelope = zeros(1,length(range_vector));
                if reset_steps
                    past_vals = zeros(1,2);
                    possible_peak=0;
                    local_peak_counter = 0;
                    forward_step_counter = 0;
                    backward_step_counter = 0;
                    min_dist = 200;
                    peak_loc = zeros(1,2); % Store two peak locations to compare them.
                end

                if Torso_val>0.95
                    counter_walking = counter_walking+1;
                    Torso_freq = f_STFT(f_ind);
                    Torso_speed = v_STFT(f_ind);
                    if (sign(walking_speed)+sign(Torso_speed))==0 % Opposite to each other
                        walking_speed=0;
                        counter_walking = 0;
                        reset_steps = 1;
                    else
                        if abs(Torso_speed)>0.1 %threshold for velocity
                            reset_steps = 0;
                            walking_speed = walking_speed+(1/counter_walking)*(Torso_speed - walking_speed);
                            ph_torso_loc.XData=range_loc_torso;
                            ph_torso_loc.YData=Torso_freq;
                            STFT_new_thresholded(:,((range_vector<range_loc_torso-0.75) | (range_vector>range_loc_torso+0.75))) =0;

                            cum_amp_distr = cumsum(STFT_new_thresholded,1)./sum(STFT_new_thresholded,1);
                            cum_amp_distr(isnan(cum_amp_distr))=0;
                            for n=1:length(range_vector)
                                freq_range = f_STFT(cum_amp_distr(:,n)>=0.025 & cum_amp_distr(:,n)<=0.975);
                                if freq_range
                                    upper_envelope(n) = freq_range(end);
                                    lower_envelope(n) = freq_range(1);
                                end
                            end
                            % Smooth the envelopes
                            upper_envelope = movmean(upper_envelope,3);
                            lower_envelope = movmean(lower_envelope,3);

                            if sign(walking_speed)==1
                                past_vals = circshift(past_vals,-1);
                                past_vals(2) = mean(maxk(upper_envelope,3)); % Current value
                                if ~possible_peak && (past_vals(2)>=past_vals(1))
                                    possible_peak = 1;
                                end
                                if possible_peak && (past_vals(2)<past_vals(1))
                                    local_peak_counter = local_peak_counter+1;
                                    if ((info_1 - peak_loc(1))<=min_dist) && (any(peak_loc))
                                        if past_vals(1)>peak_loc(2)
                                            peak_loc = [info_1,past_vals(1)];
                                        end
                                    elseif ((info_1 - peak_loc(1))>min_dist) && (any(peak_loc))
                                        forward_step_counter=forward_step_counter+1;
                                        peak_loc = [info_1,past_vals(1)];
                                    end
                                    
                                    if ~any(peak_loc)
                                        peak_loc = [info_1,past_vals(1)];
                                        forward_step_counter=forward_step_counter+1;
                                    end
                                    possible_peak=0;
                                end

                            elseif sign(walking_speed)==-1
                                past_vals = circshift(past_vals,-1);
                                past_vals(2) = mean(maxk(abs(lower_envelope),3)); % Current value
                                if ~possible_peak && (past_vals(2)>=past_vals(1))
                                    possible_peak = 1;
                                end
                                if possible_peak && (past_vals(2)<past_vals(1))
                                    local_peak_counter = local_peak_counter+1;
                                    if ((info_1 - peak_loc(1))<=min_dist) && (any(peak_loc))
                                        if past_vals(1)>peak_loc(2)
                                            peak_loc = [info_1,past_vals(1)];
                                        end
                                    elseif ((info_1 - peak_loc(1))>min_dist) && (any(peak_loc))
                                        backward_step_counter=backward_step_counter+1;
                                        peak_loc = [info_1,past_vals(1)];
                                    end

                                    if ~any(peak_loc)
                                        peak_loc = [info_1,past_vals(1)];
                                        backward_step_counter=backward_step_counter+1;
                                    end
                                    possible_peak=0;
                                end
                            end


                        else
                            walking_speed=0;
                            counter_walking=0;
                            reset_steps = 1;
                        end
                    end
                else
                    counter_walking = 0;
                    walking_speed = 0;
                    reset_steps = 1;
                end

                % Update data in handles for plotting

                % Plot received signal
                ph_signal.YData = abs(Data_Matrix_1(3,:))';
                peak_val = max([0.03 max(abs(Data_Matrix_1(3,:))) + 0.01]);
                ylim(ax_signal,[0 peak_val]);

                % Plot doppler information
%                 ph_doppler.CData = 20*log10(abs(conj(STFT))+eps);
                ph_doppler.CData = STFT_new;
                ph_upper_envelope.YData = upper_envelope;
                ph_lower_envelope.YData = lower_envelope;

                % Plot received signal history
                ph_history.CData = Data_Matrix_1_plot(1:plot_steps:end,:);
                %                 set(ax_history,'CLim',[0 peak_val+0.01]); % Set color limits of pcolor
                
                if sign(walking_speed)==1
                    general_title.String = sprintf(['{\\bfRecording Radar Data} - Date: ' datestr(now, 'dd/mm/yyyy'),' , Time: ',datestr(now,'HH:MM:SS:FFF') newline ...
                                                        'Direction: APPROACHING , Estimated walking velocity: ' sprintf('%.3f',walking_speed) ' m/s' newline ...
                                                        'Estimated Number of Forward Steps: ' num2str(forward_step_counter)]);
                elseif sign(walking_speed)==-1
                    general_title.String = sprintf(['{\\bfRecording Radar Data} - Date: ' datestr(now, 'dd/mm/yyyy'),' , Time: ',datestr(now,'HH:MM:SS:FFF') newline ...
                                                    'Direction: RECEDING    , Estimated walking velocity: ' sprintf('%.3f',walking_speed) ' m/s' newline ...
                                                    'Estimated Number of Backward Steps: ' num2str(backward_step_counter)]);
                else
                    general_title.String = sprintf(['{\\bfRecording Radar Data} - Date: ' datestr(now, 'dd/mm/yyyy'),' , Time: ',datestr(now,'HH:MM:SS:FFF') newline ...
                                                    'Direction: NO MOTION   , Estimated walking velocity: ' sprintf('%.3f',walking_speed) ' m/s' newline ...
                                                    'Estimated Number of Backward Steps: 0']);
                end
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