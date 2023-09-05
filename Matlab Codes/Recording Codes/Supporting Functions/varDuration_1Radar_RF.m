% Recording RF data for 1 Radar for a set duration of recording.

function [Data_Matrix_1,actual_fps,frame_start,frame_stop,...
         total_time,bin_length,dt,frame_axis,frame_count,range_axis,...
         frame_area_tot,range_vector_length,frame_stamp_vec,queue,time_toc,recording_actual_time]=...
            varDuration_1Radar_RF(COMPORT_1,...
                visualise_data,visualise_fps,calibration,cPLL_div,...
                PRF,iterations,pulsesperstep,dac_min,dac_max,fps,frame_area,...
                downconversion,freqRes_STFT)

% Load the library
Lib = ModuleConnector.Library;
% Display the functions available in the library - not necessary
% Lib.libfunctions;

% Create ModuleConnector object
% Log level of 0 is used below. This means that in case there is an error
% in the communication with radar we cannot continue execution after this.
% See: https://en.wikipedia.org/wiki/Syslog
mc_1 = ModuleConnector.ModuleConnector(COMPORT_1,0); %Radar 1

% Get XEP (XeThru Embedded Platform) interface
xep_1 = mc_1.get_xep(); %Radar 1

% Display systems' info
fprintf('\n-------------------------------------------------------------\n')
disp('Radar Information:')
disp(['FirmWareID = ' xep_1.get_system_info(2)]);
disp(['Version = ' xep_1.get_system_info(3)]);
disp(['Build = ' xep_1.get_system_info(4)]);
disp(['SerialNumber = ' xep_1.get_system_info(6)]);
disp(['VersionList = ' xep_1.get_system_info(7)]);

fprintf('\n')

% Clear message buffers
while xep_1.peek_message_data_float > 0
    xep_1.read_message_data_float();
end

%% Configure radar chip with x4driver through XEP interface

% First initialize chip
xep_1.x4driver_init(); % Radar 1
pause(0.1); % For stability

xep_1.x4driver_set_downconversion(downconversion);

% Maximum FPS based on radar settings above:
duty_cycle = 0.7;
FPS_max = (PRF*duty_cycle)/(iterations*pulsesperstep*(dac_max-dac_min+1));

% Note: In case fps>FPS_max, allow user to enter a valid fps value.
if fps>FPS_max
    disp('The frame rate set is higher than maximum configurable frame rate.');
    disp(['Maximum configurable frame rate is ' num2str(FPS_max) ' Hz']);
    prompt = "\nChoose a frame rate value lower than this (in Hz):   ";
    fps = input(prompt);
    while fps>FPS_max
        disp(['The frame rate set (' num2str(fps) ' Hz) is higher than maximum configurable frame rate.']);
        disp(['Maximum configurable frame rate is ' num2str(FPS_max) ' Hz']);
        prompt = "\nChoose a frame rate value lower than this (in Hz):   ";
        fps = input(prompt);
    end
end

% Set the maximum possible frame rate for both radars
xep_1.x4driver_set_prf_div(cPLL_div);
xep_1.x4driver_set_iterations(iterations);
xep_1.x4driver_set_pulsesperstep(pulsesperstep);
xep_1.x4driver_set_dac_min(dac_min);
xep_1.x4driver_set_dac_max(dac_max);

% Set frame area offset
% This is a HW (hardware) dependent constant needed to adjust for the propagation of
% the signal from the X4 chip through the antenna feed and out in the air.
% The value given here is the value found for the HW platforms in this
% example, but a different platform might have a different offset, meaning
% the actual range of objects will appear at an offset.
xep_1.x4driver_set_frame_area_offset(0.18);

% Set frame area
xep_1.x4driver_set_frame_area(frame_area(1), frame_area(2));

% Read back actual frame area
[frame_start,frame_stop] = xep_1.x4driver_get_frame_area();

% % Sampling Rate of Radar:
% fs = 23.328e9;
% 
% % Set frequency band used by the radar.
% %   frequency_band = 3 : 6.0 to 8.5 Ghz range (default mode)
% %   frequency_band = 4 : 7.2 to 10.2 GHz range
% frequency_band = 3;
% if frequency_band ==3
%     fc = 7.29e9; %Center frequency
% else
%     fc = 8.748e9; %Center frequency
% end

% xep.x4driver_set_tx_center_frequency(frequency_band);
% For some reason the following code returns a status of 1 (ERROR)
% [frequency, status] = xep.x4driver_get_tx_center_frequency();

%% Start the recording
pause(0.2)
% Start streaming data by setting FPS
% Note: This must be less than or equal to the maximum Frame Rate set by
% radar settings.
xep_1.x4driver_set_fps(fps); % Sets frame rate for frame streaming (in Hz).
disp('Recording data...')

%% The recording starts when fps>0.
actual_fps = xep_1.x4driver_get_fps();
fprintf('\nActual FPS = %f Hz\n',actual_fps);
disp(['Maximum FPS configured = ' num2str(FPS_max) ' Hz']);

%% Recording and Visualizing data (visualisation is optional)

% Most of the code in this section is to handle visualization of the data
% float message data from the module. Reading the data float message from
% the module is done with the command xep.read_message_data_float().
[~,Data_Matrix_1,bin_length,range_vector,total_time,info_1,frame_stamp_vec,queue,time_toc,recording_actual_time,frame_area_tot,range_vector_length]=rec_and_vis_data_1Radar_varDuration_RF(xep_1,frame_start,frame_stop,actual_fps,fps,visualise_data,visualise_fps,calibration,freqRes_STFT);

%% Stop radar and close connection - Copy these in Command Window in case of error 
% Use especially the xep.module_reset() function in case of error.

% Stop streaming by setting FPS = 0.
xep_1.x4driver_set_fps(0);

% Reset module
xep_1.module_reset(); % Equivalent to unplugging and replugging the radar
fprintf('Radar is reset and is ready to record new data.\n\n')
pause(0.5)

% Clean up.
clear mc_1;
clear xep_1;
clear recorder_1;
Lib.unloadlib;
clear Lib;

%Pause for 2 seconds - Leaves plot open for 2 seconds after end of
%recording.
pause(2);

%% Check Dropped frames
dt = 1/actual_fps;
frame_axis = (0:info_1-1)*dt; %Frame axis in seconds (s)

% figure;
% % Visualise which frames are dropped from Radar 1
% % 1 - Not dropped frame, 0 - Dropped frame.
% frame_stamp_Radar1 = (frame_stamp_vec(1,:)==0);
k = find(frame_stamp_vec(1,:)==0);
% subplot(1,2,1)
% stem(frame_axis,frame_stamp_Radar1,'LineStyle','none');
% title('Radar 1 - Visualising dropped frames'); xlabel('Frame index');
% ylim([-0.1 1.1]);
% yticks([0 1]); yticklabels({'Frame not Dropped','Frame Dropped'})
% xlabel('Time (s)');

fprintf('\n')
if ~isempty(k)
    disp(['Number of frames dropped from Radar: ' num2str(length(k))]);
    disp(['Percentage of frames dropped from Radar: ' num2str(100*(length(k)/size(frame_stamp_vec,2))) ' %']);
    disp(['Indices of dropped Frames: ' num2str(k)]);
else
    disp('Radar dropped 0 frames');
end
pause(0.1)

fprintf('\n')

%% Extract total number of frames and construct range axis.
range_axis = range_vector; %Range axis in meters (m)
%To find range_axis in terms of bin number:
%range_axis=(1:data_length/2);
frame_count = info_1;
end