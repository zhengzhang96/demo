function downconversion = extract_rec_file_info_empty_room(filename)

% The first part of filename is Radar_Data_ which is 11 characters
filename_cropped = filename(12:end);

% The next part indicates if data is RF or IQ.
downconv_str = filename_cropped(1:2);
downconversion = contains(downconv_str,'IQ'); % If true, data is IQ, else is RF

end