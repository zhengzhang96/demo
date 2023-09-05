function replay_SavedData_RF(Data_Matrix_1,frame_axis,range_axis,actual_fps,visualise_fps,subj_name,rec_date_ID,rec_time_ID,recording_actual_time,total_time,frame_start,frame_stop,frame_stamp_vec,freqRes_STFT)

% Compute steps in samples to visualise data (or update waitbar) at the
% fixed fps defined by visualise_fps.
visualise_steps = ceil(actual_fps/(visualise_fps)); %in samples

dt = 1/actual_fps;

% Compute dropped frames (number and percentage)
k = find(frame_stamp_vec==0);

% Display basic information about recorded data
clc
disp('------------------------------------------------');
disp(['Information about saved data: ' newline]);
disp(['Subject''s Name/ID:  ' subj_name]);
disp('Data type:  Radio-Frequency (RF)');
disp(['Recording Date (File ID):  ' rec_date_ID]);
disp(['Recording Time (File ID):  ' rec_time_ID]);
disp(['Actual Recording Start Time:  ' recording_actual_time{1}]);
disp(['Actual Recording Stop Time:  ' recording_actual_time{2}]);
disp(['Frame Rate:  ' num2str(actual_fps) ' Hz']);
disp(['Total Recording Length:  ' num2str(total_time) ' s']);
disp(['Radar Frame Area:  ' num2str(round(frame_start,3)) ' m - ' num2str(round(frame_stop,3)) ' m']);
if ~isempty(k)
    disp(['Number of frames dropped:  ' num2str(length(k)) '(' num2str(100*(length(k)/size(frame_stamp_vec,2))) ' %)']);
else
    disp('Number of frames dropped:  0 (0%)');
end
disp(['------------------------------------------------' newline]);
disp('Replaying Saved Data...');

f_0 = 7.29e9; %Center frequency
c=3e8; %Speed of light in m/s.

% Settings for RF to IQ conversion
t_fast_conversion = (2/c)*(range_axis)';
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
range_vector_IQ = range_axis(1:q:end);
Data_Matrix_1_IQ = single(NaN(size(Data_Matrix_1,1),length(range_vector_IQ)));

% Settings for Doppler Spectrum plot
if (rem(actual_fps,freqRes_STFT)/freqRes_STFT)>=0.5 || (rem(actual_fps,freqRes_STFT)/freqRes_STFT)==0
    f_num_STFT = floor(actual_fps/freqRes_STFT)+1;
else
    f_num_STFT = floor(actual_fps/freqRes_STFT);
end
win_size = 100;
window = kaiser(win_size,15);
window_mtx = repmat(window,1,size(Data_Matrix_1_IQ,2)); % Arrange window as a matrix to be applied to all range bins simultaneously.
noverlap = 75;
% determine frame shift (step) in samples
step = win_size - noverlap;
% win_centre_ind = floor(win_size/2)+1 : step : size(Data_Matrix_1,1)-(floor(win_size/2))+1; % Corresponding to the center of the window
win_upper_ind = win_size:step:size(Data_Matrix_1_IQ,1);
win_lower_ind = win_upper_ind-win_size+1;

f_STFT = (-actual_fps/2):actual_fps/f_num_STFT: (actual_fps/2) - (actual_fps/f_num_STFT);
f_0 = 7.29e9; %center frequency of carrier sinusoid.
c=3e8; %Speed of light in m/s.
lambda = c/f_0;
v_STFT = (lambda/2)*f_STFT;
% t_STFT = win_upper_ind*dt;

fh_RF = figure('Name','Replaying Saved RF signals','NumberTitle','off','Position',[218 62 1100 700]);
clf(fh_RF); %clf(fig) clears the single figure with handle fig.
fig_w = 90; %width of figure (Distance between the right and left inner edges of the figure.);
fig_h = 25; %height of figure (Distance between the top and bottom inner edges of the window.);
ButtonHandle = uicontrol('Style', 'PushButton', ...
        'String', 'Stop Playback', ...
        'Callback', 'delete(gcbf)','Position',[30 640 fig_w fig_h]);

% Plotting the received signal
subplot(4,4,[1,2,3,4])
ph_signal_RF = plot(NaN,NaN); %ph stands for plot handle
th_signal = title(''); %th stands for title handle
th_signal.String = 'Saved Signal ({\color{blue}RF},{\color{red}IQ})';
ax_signal = gca; %ax stands for axes handle
ylim(ax_signal,[-0.02 0.02]);
xlim([0.2 range_axis(end)]);
grid on;
xlabel('Range [m]');
ylabel('Amplitude (AU)')
ph_signal_RF.XData = range_axis; %x-axis values
ph_signal_RF.YData = Data_Matrix_1(1,:)';
hold on
ph_signal_IQ = plot(NaN,NaN); %ph stands for plot handle
ph_signal_IQ.XData = range_vector_IQ; %x-axis values
ph_signal_IQ.YData = Data_Matrix_1_IQ(1,:)';
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
xlabel('Range (m)'); ylabel('Elapsed Time (s)');
ph_history.XData = range_axis; %x-axis values
xlim([range_axis(1) range_axis(end)]);
ph_history.YData = frame_axis;
ylim([frame_axis(1) frame_axis(end)]);
ph_history.CData = zeros(size(Data_Matrix_1));

general_title = sgtitle ('');
general_title.FontSize = 12;

drawnow;


% Display the recorded data
for i=1:size(Data_Matrix_1,1) % run all frames one by one
    
    if ~ishandle(ButtonHandle)
        break;
    end
    
    % Convert data from IQ to RF
    RF_data_complex = (Data_Matrix_1(i,:))'.*complex_fact;
    % Apply filter - Note: filtfilt operates along the first array dimension with size greater than 1.
    RF_data_complex_filt = FiltFiltM(b_lp,a_lp,RF_data_complex);
    IQ_data = RF_data_complex_filt(1:q:end,:)';

    Data_Matrix_1_IQ(i,:) = conj(IQ_data);

    % Apply STFT in real-time - FFT on data*window function.
    [TF,loc]=ismember(i,win_upper_ind);
    if TF
        % Multiply segment of signal by window function
        Data_Matrix_1_window = Data_Matrix_1_IQ(win_lower_ind(loc):win_upper_ind(loc),:).*window_mtx;
        % Compute the Fourier Transform of this segment of data
        STFT = fftshift(fft(Data_Matrix_1_window,f_num_STFT,1),1); % fft along first dimension.
    end
    
    if ~mod(i,visualise_steps) && i>=win_upper_ind(1) && i<=win_upper_ind(end)% Plot at the defined visualisation frame rate

        % Plot received signal
        ph_signal_RF.YData = Data_Matrix_1(i,:)';
        ph_signal_IQ.YData = abs(Data_Matrix_1_IQ(i,:))';
        peak_val = max(Data_Matrix_1(i,:));
        ylim(ax_signal,[min([-0.02 (-peak_val - 0.01)]) max([0.02 (peak_val + 0.01)])]);

        % Plot doppler information
        ph_doppler.CData = 20*log10(abs(conj(STFT))+eps);

        % Plot received signal history
        Data_Matrix_1_plot = [Data_Matrix_1(1:i,:);zeros(size(Data_Matrix_1,1)-i,size(Data_Matrix_1,2))];
        ph_history.CData = abs(Data_Matrix_1_plot);
        set(ax_history,'CLim',[0 max([0.02 max(max(abs(gather(Data_Matrix_1(1:i,:))))) + 0.005])]); % Set color limits of pcolor

        general_title.String = sprintf('Saved Data - Elapsed Time = %.2f s',round(frame_axis(i)+dt,2));
        
        pause(1/visualise_fps)
        drawnow;
    end

end
end