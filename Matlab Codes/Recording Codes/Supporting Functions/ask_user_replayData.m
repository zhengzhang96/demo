function replay = ask_user_replayData()

confirm = 0;

while ~confirm

    % Ask user if they want to replay data
    prompt = '\nWould you like to replay recorded data (1-Yes, 0-No)?   ';
    replay = input(prompt);
    while replay~=0 && replay~=1
        disp('Invalid Input');
        prompt = '\nWould you like to replay recorded data (1-Yes, 0-No)?   ';
        replay = input(prompt);
    end

    % Display settings to user for confirmation
    if replay
        disp('You have selected to replay the recorded data.');
    else
        disp('The recorded data will not be replayed.');
    end

    prompt = "\nConfirm (1-Yes, 0-No)?   ";
    confirm = input(prompt);
    while confirm~=0 && confirm~=1
        disp('Invalid Input');
        prompt = "\nConfirm (1-Yes, 0-No)?   ";
        confirm = input(prompt);
    end

end

end