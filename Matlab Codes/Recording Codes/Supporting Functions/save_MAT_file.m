function save_MAT_file(filename,recording_actual_time,Data_Matrix_1,frame_stamp_vec,frame_axis,range_axis,calibration,actual_fps,visualise_fps,frame_start,frame_stop,total_time,bin_length,dt,frame_count,freqRes_STFT)

%Storing data as .mat file
if ~exist('./Recordings MAT files', 'dir')
    %Save files in a folder (folder created if non existent).
    mkdir('./Recordings MAT files')
end

% Complete the file name for .mat files
filename_MAT = append('./Recordings MAT files/Radar_Data_',filename,'.mat');

% Attach actual recording start and stop times for future reference
recording_start_time=recording_actual_time{1};
recording_stop_time=recording_actual_time{2};

disp('Saving Data as a .MAT file ...')

% Convert values of data matrices to 32-bit
Data_Matrix_1 = single(Data_Matrix_1);
frame_stamp_vec = single(frame_stamp_vec);
frame_axis = single(frame_axis);
range_axis = single(range_axis);

% Saving WITH compression makes saving process very slow. To remove
% compression, use '-nocompression'. However, files are much larger now
% in size. (We can do compression after trials if needed).
% Note: If a variable is >=2GB, we have an error/warning from MATLAB.
% In this case there are two options:
% 1) use version '-v7.3' (slower than default)
% 2) split file into multiple files.

if ~calibration
    save(fullfile(filename_MAT),'Data_Matrix_1','actual_fps','visualise_fps',...
        'frame_start','frame_stop','total_time','bin_length','dt','frame_axis','freqRes_STFT',...
        'frame_count','frame_stamp_vec','range_axis','recording_start_time','recording_stop_time','-nocompression')
else
    save(fullfile(filename_MAT),'Data_Matrix_1','frame_axis','range_axis','-nocompression')
end

disp('Saving Data COMPLETE')

end