function save_empty_room_file(filename,Data_Matrix_1,frame_axis,range_axis)

%Storing data as .mat file
if ~exist('./Empty Room Recordings MAT files', 'dir')
    %Save files in a folder (folder created if non existent).
    mkdir('./Empty Room Recordings MAT files')
end

% Complete the file name for .mat files
filename_MAT = append('./Empty Room Recordings MAT files/Radar_Data_',filename,'.mat');

disp('Saving Empty Room Data as a .MAT file ...')

% Convert values of data matrices to 32-bit
Data_Matrix_1 = single(Data_Matrix_1);
frame_axis = single(frame_axis);
range_axis = single(range_axis);

save(fullfile(filename_MAT),'Data_Matrix_1','frame_axis','range_axis','-nocompression')


disp('Saving Empty Room Data COMPLETE')

end