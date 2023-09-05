function [Data_Matrix_1,actual_fps,frame_start,...
          frame_stop,total_time,bin_length,dt,frame_axis,...
          frame_count,range_axis,frame_area_tot,range_vector_length,frame_stamp_vec,...
          queue,time_toc,recording_actual_time,calibration,filename,subj_name,time] = ...
               select_recording_fun(...
                                   set_duration,downconversion,COMPORT_1,...
                                   visualise_fps,cPLL_div,PRF,iterations,pulsesperstep,dac_min,dac_max,fps,...
                                   frame_area,freqRes_STFT)

% Ask if user wants to record empty room data for calibration purposes
calibration = ask_user_calibration();

if set_duration % For set duration of the recording
   
    if ~calibration
        visualise_data=1;
        total_duration = ask_user_info_setDuration();
        subj_name = GenerateSubjectName();
        % Generate file name (for Vicon and .MAT files):
        time = datestr(now, 'yyyymmdd_HHMMSS'); %Adding Time mark on file
        if downconversion
            filename = sprintf('IQ_%s',time);
        else
            filename = sprintf('RF_%s',time);
        end
        %Adding subject's name on file.
        filename = append(filename,'_',subj_name);
    else
        visualise_data=0;
        total_duration=20; %20 s for empty room recording.
        if downconversion
            filename = sprintf('IQ_Calibration');
        else
            filename = sprintf('RF_Calibration');
        end
    end

    if downconversion % For Recording Baseband (IQ) data
        [Data_Matrix_1,actual_fps,frame_start,...
            frame_stop,total_time,bin_length,dt,frame_axis,...
            frame_count,range_axis,frame_area_tot,range_vector_length,frame_stamp_vec,...
            queue,time_toc,recording_actual_time]=...
            setDuration_1Radar_IQ(COMPORT_1,...
            visualise_data,visualise_fps,calibration,cPLL_div,...
            PRF,iterations,pulsesperstep,dac_min,dac_max,fps,frame_area,...
            downconversion,total_duration,freqRes_STFT);
    else % For Recording RF data
        [Data_Matrix_1,actual_fps,frame_start,frame_stop,...
            total_time,bin_length,dt,frame_axis,frame_count,range_axis,...
            frame_area_tot,range_vector_length,frame_stamp_vec,queue,time_toc,recording_actual_time]=...
            setDuration_1Radar_RF(COMPORT_1,...
            visualise_data,visualise_fps,calibration,cPLL_div,...
            PRF,iterations,pulsesperstep,dac_min,dac_max,fps,frame_area,...
            downconversion,total_duration,freqRes_STFT);
    end

else % For variable duration of the recording
    if ~calibration
        visualise_data=1;
        subj_name = GenerateSubjectName();
        % Generate file name (for Vicon and .MAT files):
        time = datestr(now, 'yyyymmdd_HHMMSS'); %Adding Time mark on file
        if downconversion
            filename = sprintf('IQ_%s',time);
        else
            filename = sprintf('RF_%s',time);
        end
        %Adding subject's name on file.
        filename = append(filename,'_',subj_name);
    else
        visualise_data=0;
        total_duration=20; %20 s for empty room recording.
        if downconversion
            filename = sprintf('IQ_Calibration');
        else
            filename = sprintf('RF_Calibration');
        end
    end

    if downconversion % For Recording Baseband (IQ) data
        if ~calibration
            [Data_Matrix_1,actual_fps,frame_start,frame_stop,...
                total_time,bin_length,dt,frame_axis,frame_count,range_axis,...
                frame_area_tot,range_vector_length,frame_stamp_vec,queue,time_toc,recording_actual_time]=...
                varDuration_1Radar_IQ(COMPORT_1,...
                visualise_data,visualise_fps,calibration,cPLL_div,...
                PRF,iterations,pulsesperstep,dac_min,dac_max,fps,frame_area,...
                downconversion,freqRes_STFT);
        else
            [Data_Matrix_1,actual_fps,frame_start,...
                frame_stop,total_time,bin_length,dt,frame_axis,...
                frame_count,range_axis,frame_area_tot,range_vector_length,frame_stamp_vec,...
                queue,time_toc,recording_actual_time]=...
                setDuration_1Radar_IQ(COMPORT_1,...
                visualise_data,visualise_fps,calibration,cPLL_div,...
                PRF,iterations,pulsesperstep,dac_min,dac_max,fps,frame_area,...
                downconversion,total_duration,freqRes_STFT);
        end
    else % For Recording RF data
        if ~calibration
            [Data_Matrix_1,actual_fps,frame_start,frame_stop,...
                total_time,bin_length,dt,frame_axis,frame_count,range_axis,...
                frame_area_tot,range_vector_length,frame_stamp_vec,queue,time_toc,recording_actual_time]=...
                varDuration_1Radar_RF(COMPORT_1,...
                visualise_data,visualise_fps,calibration,cPLL_div,...
                PRF,iterations,pulsesperstep,dac_min,dac_max,fps,frame_area,...
                downconversion,freqRes_STFT);
        else
            [Data_Matrix_1,actual_fps,frame_start,...
                frame_stop,total_time,bin_length,dt,frame_axis,...
                frame_count,range_axis,frame_area_tot,range_vector_length,frame_stamp_vec,...
                queue,time_toc,recording_actual_time]=...
                setDuration_1Radar_RF(COMPORT_1,...
                visualise_data,visualise_fps,calibration,cPLL_div,...
                PRF,iterations,pulsesperstep,dac_min,dac_max,fps,frame_area,...
                downconversion,total_duration,freqRes_STFT);
        end
    end
end

% Note: we use different filename for Vicon that when we save .mat files.
% This is because the .mat filename includes folder path.
% Do not use the following special characters in a Vicon trial name:
% \ . / , < ? > * : " | % $ 

end