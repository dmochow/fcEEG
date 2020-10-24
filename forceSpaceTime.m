function dataOut = forceSpaceTime(dataIn)
% im assuming this is what Loq D does in his function
if size(dataIn,1) < size(dataIn,2)
    % all good, channels(1) should be less than samples(2)
    dataOut = dataIn;
elseif size(dataIn,1) > size(dataIn,2)
    % not good, channels(1) should be less than samples(2)
    dataOut = permute(dataIn,[2 1 3]);
end