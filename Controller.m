%clear all;
%clc;

% initiate Glucose Model
GlucoseModelInit;

% simulation time (24 hrs)
tf = 24.*60/15;

%state and action definitions
%s= glucose;
%a= insulin;

% Increase insulin by Standard Treatment Policy for injections
% insulinBolus=(((glucose-Ctarget)/CF)+(carbs/CR))

%glucose= current
%Ctarget= Target Glucose
%CF= Correction Factor
%carbs= meal
%CR= Carb Ratio

%q table
% s = 70:10:200;
s = [ 50 89 110 131 180 ]; %  hypoglycaemia, normoglycaemia, hyperglyceamia
a = round( 1000*(0:0.2:1) );
% Q = abs(randn( length(s), length(a) ));
Q = 0.1 + 0.6.*rand(length(s), length(a) );
% Q(s==150,a==1000) = 1;
% Q(s==60,a==0) = 1;
% Q=Q_Good;

% %Cumulative Q
% Qcm= sum(Q(:))

%start Q-Learning

% %random intial glucose
% BG_0= round(140*rand());

% % match glucose to state
% diff = abs( BG_0 - s );
% [dum,id] = min( diff );
% s_current = s( id );

% match insulin to action [not sure if necessary]
%diff1 = abs( y - a );
%[dum1,id1] = min( diff1 );
%a_current = a( id1 );

%having some problems with this
%y=RLpolicy(BG_0);

mdl = 'GlucoseModelSubsystem2016b';


%Update Q-Table
max_iterations = 24*60/15; %15 minute intervals or not depending on model rates
alpha= 0.33;
gamma= 0.99;
epsilon = 0.75;

mdl = 'GlucoseModelSubsystem2016b';
open_system( mdl );

% initial step





for k = 1 : 100
    Q0 = Q;
    TIME = [];
    BGL = [];
    INS = [];
    
    tstop = 1;
    set_param( mdl, ...
        'SaveFinalState', 'on', ...
        'FinalStateName', [ mdl 'SimState' ], ...
        'SaveCompleteFinalSimState', 'on' );
    
    simOut = sim( mdl, 'StopTime', num2str(tstop) );
    
    
    for i = 1 : max_iterations
        
        % simulation step
        
        InitState = simOut.get( [ mdl 'SimState' ] );
        tstop = tstop + 15;
        set_param( mdl, ...
            'LoadInitialState', 'on', ...
            'InitialState', 'InitState' );
        simOut = sim( mdl, 'StopTime', num2str(tstop) );
        % this gets the new state (glucose level) from simulation
        tmp = simOut.get('glucose');
        BGL0 = 18.*tmp(tstop-15);
        BGL(i,:) = 18.*tmp(tstop);
        
        tmp1= simOut.get('insulin');
        INS(i,:) = tmp1(tstop);
        
        tmp3 = simOut.get('tout');
        TIME(i,:) = tmp3(tstop);
        
        %id_s - index of the current state
        
        %match new glucose to s
        diff_next= abs( BGL0 - s );
        [dum,id_s] = min( diff_next );
        s_now = s( id_s );
        
        diff_next= abs( BGL(i,:) - s );
        [dum,id_s_next] = min( diff_next );
        s_next = s( id_s_next );
        
%         %id_a - index of the current action
        diff_next= abs( INS(i,:) - a );
        [dum,id_a] = min( diff_next );
%         a_now = a( id_a );
%         % perform
        
        r = reward( s_next );
        %update the Q table
        Q(id_s,id_a) = Q(id_s,id_a) + alpha * ( r + gamma*max(Q(id_s_next,:)) - Q(id_s,id_a) );
        %     Q(id_s, id_a) = R(id_s, id_a) + Gamma * max(Q(id_next, a))
        if length(TIME) > 2
            hf = figure(1);
            set(hf,'Position',[ 50 50 1000 800 ]);
            subplot(211);
            plot( TIME, BGL, 'm' );
            hold on;
            line( [ TIME(1) TIME(end) ],[90 90],'LineStyle',':','Color','k');
            line( [ TIME(1) TIME(end) ],[130 130],'LineStyle',':','Color','k');
            line( [ TIME(1) TIME(end) ],[70 70],'LineStyle',':','Color','r');
            xlim( [ TIME(1) TIME(end) ] );
            ylim( [ 0 300 ] );
            
            subplot(212);
            stairs( TIME, INS, 'b' );
            xlim( [ TIME(1) TIME(end) ] );
            ylim( [ 0 1.1.*max(a) ] );

            drawnow;
        end
        
                 figure(2);
                surf(Q);
                shading interp;
                drawnow;
    end
    Q1 = Q;
    Qnorm(k) = norm( Q1 - Q0 )
    Qsum(k) = sum(Q1(:))
    set_param( mdl, ...
        'SaveFinalState', 'off', ...
        'LoadInitialState', 'off' );
    epsilon = 0.9.*epsilon;
    alpha = 0.99.*alpha;
    subplot(211); hold off;
    
end


%reward
