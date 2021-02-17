function SentMessages = SendMessagesMatrix_v1(SendMessages,p,r)
% SendMessages: a matrix filled with entries equal to 1 if it was sent last time and 0 if not 
% p: probability of losing the i-th message if the previous one was send
% r: probability of starting to send the i-th message if the previous one
%was lost
% GILBERT-ELLIOT MODEL (loss bursts)

SentMessages = SendMessages;
[size1,size2] = size(SendMessages);
for i=1:size1
    for j = 1:size2
        if(SendMessages(i,j)==1 && (rand(1) <= p) )
                SentMessages(i,j) = 0;
        elseif(SendMessages(i,j)==0 && (rand(1) <= r) )
                     SentMessages(i,j) = 1;
        end
    end
end

