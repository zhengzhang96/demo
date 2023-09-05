function total_duration = ask_user_info_setDuration()

confirm = 0;

while ~confirm

    %Total Duration of recording
    prompt = '\nTotal Duration of Recording in seconds?   ';
    total_duration = input(prompt);
    while total_duration<=0
        disp('Invalid Input');
        prompt = '\nTotal Duration of Recording in seconds?   ';
        total_duration = input(prompt);
    end

    % Display settings to user for confirmation
    clc
    disp(['Selected Duration of recordings = ' num2str(total_duration) ' s']);

    prompt = "\nConfirm settings (1-Yes, 0-No)?   ";
    confirm = input(prompt);
    while confirm~=0 && confirm~=1
        disp('Invalid Input');
        prompt = "\nConfirm settings (1-Yes, 0-No)?   ";
        confirm = input(prompt);
    end

    if ~confirm
        clc
    end

end

end