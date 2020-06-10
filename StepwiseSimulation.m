clear all 
clc

GlucoseModelInit;

%q table
s = 50:10:150;
a = 60.*(0:5:40);
Q = rand( length(s), length(a) );


mdl = 'GlucoseModelSubsystem2016b';
open_system( mdl );
set_param( mdl, ...
    'SaveFinalState', 'on', ...
    'FinalStateName', [ mdl 'SimState' ], ...
    'SaveCompleteFinalSimState', 'on' );


tstop = 1;
simOut = sim( mdl, 'StopTime', num2str(tstop) );
glucose_out(1) = simOut.glucose(2)

% this part of the code performs state transition
InitState = simOut.get( [ mdl 'SimState' ] );
tstop = tstop + 1;
set_param( mdl, ...
    'LoadInitialState', 'on', ...
    'InitialState', 'InitState' );
simOut = sim( mdl, 'StopTime', num2str(tstop) );
% this gets the new state (glucose level) from simulation
glucose_out(2) = simOut.glucose(2)







set_param( mdl, 'LoadInitialState', 'off' );