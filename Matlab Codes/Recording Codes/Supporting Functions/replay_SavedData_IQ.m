function replay_SavedData_IQ(Data_Matrix_1,frame_axis,range_axis,actual_fps,visualise_fps,subj_name,rec_date_ID,rec_time_ID,recording_actual_time,total_time,frame_start,frame_stop,frame_stamp_vec,freqRes_STFT)

% Compute steps in samples to visualise data (or update waitbar) at the
% fixed fps defined by visualise_fps.
visualise_steps = ceil(actual_fps/(visualise_fps)); %in samples

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

dt = 1/actual_fps;

% Compute dropped frames (number and percentage)
k = find(frame_stamp_vec==0);

% Display basic information about recorded data
clc
disp('------------------------------------------------');
disp(['Information about saved data: ' newline]);
disp(['Subject''s Name/ID:  ' subj_name]);
disp('Data type:  Baseband (I/Q)');
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


% Initialise figure
fh_baseband = figure('Name','Replaying Saved Baseband signals','NumberTitle','off','Position', [218 62 1100 700]);
clf(fh_baseband);
fig_w = 90; %width of figure (Distance between the right and left inner edges of the figure.);
fig_h = 25; %height of figure (Distance between the top and bottom inner edges of the window.);
ButtonHandle = uicontrol('Style', 'PushButton', ...
        'String', 'Stop Playback', ...
        'Callback', 'delete(gcbf)','Position',[30 640 fig_w fig_h]);
% Plotting the received signal
subplot(4,4,[1,2,3,4])
ph_signal = plot(NaN,NaN); %ph stands for plot handle
th_signal = title(''); %th stands for title handle
th_signal.String = 'Saved Signal';
ax_signal = gca; %ax stands for axes handle
ylim(ax_signal,[0 0.05]);
xlim([range_axis(1) range_axis(end)]);
grid on;
xlabel('Range [m]');
ylabel('Magnitude (AU)')
ph_signal.XData = range_axis; %x-axis values
ph_signal.YData = abs(Data_Matrix_1(1,:))';

% Plotting the doppler information
subplot(4,4,[7,8,11,12,15,16])
yyaxis left % Showing y-axis in Hz
ph_doppler = imagesc(0,0,0);
th_doppler = title('');
th_doppler.String = 'Range-Doppler Information';
ax_doppler = gca;
set(ax_doppler,'YDir','normal')
% set(ax_doppler,'CLim',[0.75 1]);
colormap(jet)
xlabel('Range [m]'); ylabel('Doppler Frequency (Hz)');
xlim([range_axis(1) range_axis(end)]);
ylim([f_STFT(1) f_STFT(end)]);
ph_doppler.XData = range_axis; %x-axis values
ph_doppler.YData = f_STFT; %y-axis values
ph_doppler.CData = 20*log10(eps*ones(length(f_STFT),length(range_axis))); %y-axis values
yyaxis right % Showing y-axis in m/s
ph_doppler_right = imagesc(0,0,0);
ax_doppler_right = gca;
set(ax_doppler_right,'YDir','normal')
ph_doppler_right.YData = v_STFT; %y-axis values
ph_doppler_right.XData = range_axis; %x-axis values
xlim([range_axis(1) range_axis(end)]);
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
xlabel('Range (m)'); ylabel('Elapsed Time (s)');
ph_history.XData = range_axis; %x-axis values
xlim([range_axis(1) range_axis(end)]);
ph_history.YData = frame_axis;
ylim([frame_axis(1) frame_axis(end)]);
ph_history.CData = zeros(size(Data_Matrix_1));

general_title = sgtitle ('');
general_title.FontSize = 12;

r=1;

% Display the recorded data
for i=1:size(Data_Matrix_1,1) % run all frames one by one
    
    if ~ishandle(ButtonHandle)
        break;
    end

    % Apply STFT - FFT on data*window function.
    [TF,loc]=ismember(i,win_upper_ind);
    if TF
        % Multiply segment of signal by window function
        Data_Matrix_1_window = Data_Matrix_1(win_lower_ind(loc):win_upper_ind(loc),:).*window_mtx;
        % Compute the Fourier Transform of this segment of data
        STFT = fftshift(fft(Data_Matrix_1_window,f_num_STFT,1),1); % fft along first dimension.
%         STFT(abs(STFT)<0.005) = 0;
% 
%         max(max(abs(STFT)))
        mu = mean(abs(conj(STFT)),'all');
        STFT_new = (abs(conj(STFT)).^r)./((abs(conj(STFT)).^r) + (mu.^r));
    end

    if ~mod(i,visualise_steps) && i>=win_upper_ind(1) && i<=win_upper_ind(end)
%         clf %clf(fig) clears the single figure with handle fig.

        % Plot received signal
        ph_signal.YData = abs(Data_Matrix_1(i,:))';
        peak_val = max([0.03 max(abs(Data_Matrix_1(i,:))) + 0.01]);
        ylim(ax_signal,[0 peak_val]);
% 
%         % Plot received signal history
%         Data_Matrix_1_plot = [Data_Matrix_1(1:i,:);zeros(size(Data_Matrix_1,1)-i,size(Data_Matrix_1,2))];
%         ph_history.CData = abs(Data_Matrix_1_plot);
%         set(ax_history,'CLim',[0 max([0.02 max(max(abs(gather(Data_Matrix_1(1:i,:))))) + 0.005])]); % Set color limits of pcolor

    subplot(4,4,[5,6,9,10,13,14])
    [counts,centers] = hist(20*log10(abs(conj(STFT(:)))+eps),100);
    stem(centers,log10(counts))

        % Plot doppler information
        ph_doppler.CData = abs(conj(STFT));
%         ph_doppler.CData = 20*log10(abs(conj(STFT))+eps);
%         ph_doppler.CData = STFT_new;
%         lower_lim_STFT = mean(STFT_new,'all') + 3*std(STFT_new,[],'all');
        
        general_title.String = sprintf('Saved Data - Elapsed Time = %.2f s',round(frame_axis(i)+dt,2));

        pause(1/visualise_fps)
        drawnow;
    end

end
end