% define states :
%     glucose
% define action :
%     insulin dose
% define reward :
%      1 if glucose 90-130
%      0 otherwise


% example reward function

function r = reward( s )
if (s >= 90) & (s <= 130)
    r = 1;
else
    r = 0;
end


% example q table
% q(s,a)
% s = [ 50 60 70 80 90 100 110 120 130 140 150 ]
% a = 0:5:40;
% Q = rand( length(s), length(a) );

