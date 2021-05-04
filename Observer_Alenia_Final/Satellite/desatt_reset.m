% set initial desired attitude
function desatt_reset(time_instant)
    
    global DynOpt params
    
    temp_att = zeros(3,1);
    if mod(time_instant,DynOpt.Nts) == 0
        for i=1:1
            t = 10*i*DynOpt.time(time_instant);
            temp_att = temp_att + pi/8*([sin(t); cos(t); -sin(0.7*t)]);
        end
        params.DesiredAttitude =  temp_att;
    else
        for i=1:1
            t = 10*i*DynOpt.time(time_instant-mod(time_instant,DynOpt.Nts));
            temp_att = temp_att + pi/8*([sin(t); cos(t); -sin(0.7*t)]);
        end
        params.DesiredAttitude =  temp_att;
    end
end