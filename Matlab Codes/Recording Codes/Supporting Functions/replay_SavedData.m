function replay_SavedData(Data_Matrix_1,frame_axis,range_axis,actual_fps,visualise_fps,subj_name,downconversion,rec_date_ID,rec_time_ID,recording_actual_time,total_time,frame_start,frame_stop,frame_stamp_vec,freqRes_STFT)
    % Select if recorded data is originally RF or IQ (baseband)
    if downconversion
        replay_SavedData_IQ(Data_Matrix_1,frame_axis,range_axis,actual_fps,visualise_fps,subj_name,rec_date_ID,rec_time_ID,recording_actual_time,total_time,frame_start,frame_stop,frame_stamp_vec,freqRes_STFT)
    else
        replay_SavedData_RF(Data_Matrix_1,frame_axis,range_axis,actual_fps,visualise_fps,subj_name,rec_date_ID,rec_time_ID,recording_actual_time,total_time,frame_start,frame_stop,frame_stamp_vec,freqRes_STFT)
    end
end