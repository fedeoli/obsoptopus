function ObserverUpdateJumpMap_V1_1(time_index)

global ObserverTest Agent

%updating data not linked with estimation
for k=1:ObserverTest.Nagents,
    Agent(k).xHatUKF(1:6,time_index) = Agent(k).xHatUKF(1:6,time_index) + ObserverTest.Delta_iner_ECI(1:6,time_index);
end
