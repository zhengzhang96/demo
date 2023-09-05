function rec_file_info = extract_rec_file_info(filename)

% The first part of filename is Radar_Data_ which is 11 characters
filename_cropped = filename(12:end);

% The next part indicates if data is RF or IQ.
downconv_str = filename_cropped(1:2);
downconversion = contains(downconv_str,'IQ'); % If true, data is IQ, else is RF

% Extract date and time of experimental file.
filename_cropped = filename_cropped(4:end);
date = filename_cropped(1:8); %Date
time = filename_cropped(10:15); %Time
% Convert date into DD/MM/YY
rec_Date = append(date(7:8),'/',date(5:6),'/',date(1:4));
% Convert time into HH:MM:SS
rec_Time = append(time(1:2),':',time(3:4),':',time(5:6));

% Extract subject's Name/ID:
filename_cropped = filename_cropped(17:end);
end_ind_subj = strfind(filename_cropped,'.mat')-1;
subj_id = filename_cropped(1:end_ind_subj);

% Store file information in a structure
rec_file_info.data_type = downconv_str;
rec_file_info.downconversion = downconversion;
rec_file_info.recording_date = rec_Date; 
rec_file_info.recording_time = rec_Time;
rec_file_info.subject_ID = subj_id; 

end