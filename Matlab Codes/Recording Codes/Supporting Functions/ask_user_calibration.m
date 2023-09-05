function calibration = ask_user_calibration()

confirm = 0;

while ~confirm

    % Ask if user wants to record empty room data for calibration
    prompt = '\nRecord Empty Room data for calibration (1-Yes, 0-No)?   ';
    calibration = input(prompt);
    while calibration~=0 && calibration~=1
        disp('Invalid Input');
        prompt = '\nRecord Empty Room data for calibration (1-Yes, 0-No)?   ';
        calibration = input(prompt);
    end
    
    % Display settings to user for confirmation
    clc
    if calibration
        disp('You have selected to record empty room data.');
    else
        disp('No empty room data will be recorded.');
    end
    
    prompt = "\nConfirm (1-Yes, 0-No)?   ";
    confirm = input(prompt);
    while confirm~=0 && confirm~=1
        disp('Invalid Input');
        prompt = "\nConfirm (1-Yes, 0-No)?   ";
        confirm = input(prompt);
    end
    
    if ~confirm
        clc
    end

end

end