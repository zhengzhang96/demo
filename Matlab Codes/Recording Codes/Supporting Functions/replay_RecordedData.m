function replay_RecordedData(Data_Matrix_1,frame_axis,range_axis,actual_fps,visualise_fps,subj_name,downconversion,time,recording_actual_time,total_time,frame_start,frame_stop,frame_stamp_vec,freqRes_STFT)
    % Select if recorded data is originally RF or IQ (baseband)
    if downconversion
        replay_RecordedData_IQ(Data_Matrix_1,frame_axis,range_axis,actual_fps,visualise_fps,subj_name,time,recording_actual_time,total_time,frame_start,frame_stop,frame_stamp_vec,freqRes_STFT)
    else
        replay_RecordedData_RF(Data_Matrix_1,frame_axis,range_axis,actual_fps,visualise_fps,subj_name,time,recording_actual_time,total_time,frame_start,frame_stop,frame_stamp_vec,freqRes_STFT)
    end
end