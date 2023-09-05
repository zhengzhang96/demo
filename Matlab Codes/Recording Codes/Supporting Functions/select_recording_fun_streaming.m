function select_recording_fun_streaming(COMPORT_1)

% Set downconversion - 0 = no downconversion (RF data) and 1 = downconversion (Baseband data).
downconversion = 1;

% Set visualisation frame rate (as small as possible) when diplaying data
% This only applies for non-calibration data.
visualise_fps = 10; %in Hz. (Either 5Hz, 10Hz or 20Hz).

% Frame Area of recording
frame_area(1) = 0.5; % Frame start
frame_area(2) = 6; % Frame end

% Set the frequency resolution for Doppler Spectrum.
freqRes_STFT = 1; %Frequency Resolution in Hz

% Set the frame rate in Hz - Maximum with current settings is 546.4253 Hz
fps = 500;

% X4 chip radar settings - These define the maximum possible frame rate
cPLL_div = 16; % default is 16 (in general it is set as an integer, e.g. 15.5 -> 16).
               % The valid divider range is 6-255.
PRF = (243e6)/cPLL_div; % PRF is a fraction of Common PLL clock which is 243MHz.
iterations = 8; % Valid values are 1-255 and should be divisible by 4.
pulsesperstep = 16; % Valid values are 1-65535.
dac_min = 949; % 11-bit DAC - DAC range from 0 to 2047
dac_max = 1100;

% Ask if user wants to record empty room data for calibration purposes
% (important before start streaming data)
calibration = ask_user_calibration();

% Considering only variable duration of the recording (more useful for streaming)
if ~calibration
    visualise_data=1;
else
    visualise_data=0;
    total_duration=60; %60 s (1 minute) for empty room recording.
    filename = sprintf('IQ_Calibration');
end

% Recording Baseband (IQ) data
    if ~calibration
        varDuration_1Radar_IQ_stream(COMPORT_1,...
            visualise_data,visualise_fps,calibration,cPLL_div,...
            PRF,iterations,pulsesperstep,dac_min,dac_max,fps,frame_area,...
            downconversion,freqRes_STFT);
    else
        [Data_Matrix_1,~,~,~,~,~,~,frame_axis,~,range_axis,~,~,~,~,~,~]=...
            setDuration_1Radar_IQ(COMPORT_1,visualise_data,visualise_fps,...
            calibration,cPLL_div,PRF,iterations,pulsesperstep,dac_min,...
            dac_max,fps,frame_area,downconversion,total_duration,freqRes_STFT);
    end

if calibration
    % Save empty room data in a .mat file - useful for empty room subtraction
    save_empty_room_file(filename,Data_Matrix_1,frame_axis,range_axis)
end

end