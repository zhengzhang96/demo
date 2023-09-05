function subj_name = GenerateSubjectName()

persistent subj_name_counter

if isempty(subj_name_counter)
    subj_name_counter = 1;
else
    subj_name_counter = subj_name_counter + 1;
end

subj_name = sprintf('DemoRecording%d',subj_name_counter);
end
