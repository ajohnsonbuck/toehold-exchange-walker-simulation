clear all;
close all;

% Performs kinetic Monte Carlo simulation of single-molecule Forster resonance
% energy transfer (smFRET) measurements of a DNA toehold-exchage walker as
% described in Li, Johnson-Buck, et al. Nature Nanotechnology 13, p. 723–729 (2018)

% The model models toehold binding and dissociation as well as branch migration
% processes using a kinetic Monte Carlo approach.

% The output is a FRET-versus-time trajectory for a single simulated walker
% over a user-specified measurement interval.

%% User-defined parameters-----------------------------------

nruns = 1;
nfootholds = 3; % Number of footholds in system
a = 8; % Size of toehold, in nucleotides
b = 13; % Size of branch migration domain, in nucleotides
tfinal = 5; % Final time, in seconds
tsegmentlength = 0.5; % Segment size to break simulation into
itime = 0.1; % Camera/measurement integration time, in s

%% Main body of code------------------------------------------

if tsegmentlength < tfinal
    tsegments = floor(tfinal/tsegmentlength);
else
    tsegmentlength = tfinal;
end

plot_average = 1;
calculate_FRET = 0;
%for 5 bp
FRETstates = [0.994680455	0.999271531	0.99725995	0.990363291	0.972421301	0.935129155	0.87112115	0.778780151	0.665394826	0.544901017	0.431383206	0.333818064	0.255142884	0.19420019	0.14804462	0.081165667	0.081165667	0.081165667	0.081165667	0.081165667	0.081165667	0.081165667	0.081165667	0.081165667	0.081165667	0.081165667	0.081165667	0.081165667	0.081165667	0.081165667	0.003662809];

kbind = 300;
kbindratio = 1; % Ratio between binding to F1 and F1' footholds
kdissoc = 3E6*exp(-2.031*a);
% kdissoc = 20; % 5 bp toehold, empirical
% kdissoc = 17.5; % 6 bp toehold, empirical
% kdissoc = 13.3; % 7 bp toehold, empirical
% kdissoc = 0.8889; % 8 bp toehold, empirical
% kbm = 1000;
% kbm = 3000; % Rate constant for branch migration
kbm = 10000;
% kbmi = 400; % Rate constant for initiating branch migration
kbmi = 1400;

nstates = nfootholds + (nfootholds-1)*(b+1);

% Create rate constant matrix, initializing such that all states have rate
% constants kbm
kmat = zeros(nstates,2)+kbm;

% Set binding rate constants
for m = 1:(b+2):nstates
    kmat(m,1) = kbind;
    kmat(m,2) = kbind;
end

kmat(b+3,2) = kbind/kbindratio;
kmat(1,1) = kbind/kbindratio;

% Set reverse rate of first state and forward rate of last state == 0
kmat(1,2) = 0;
kmat(end,1) = 0;

% Set rates for toehold dissociation steps and branch migration initiation
for m = 2:(b+2):nstates
    kmat(m,2) = kdissoc;
    kmat(m,1) = kbmi;
end

for m = b+2:(b+2):nstates
    kmat(m,1) = kdissoc;
    kmat(m,2) = kbmi;
end

steptime = 1./(kmat(:,1)+kmat(:,2));
Pfwd = kmat(:,1)./(kmat(:,1) + kmat(:,2));

for run = 1:nruns

    if mod(run,2)==0
        disp('Starting run:');
        disp(num2str(run));
    end

    m = unidrnd(2*b+5); % Determine starting state
    t = 0;
    statelistbin_all = zeros(0,1);
    binlist_all = zeros(0,1);

    tlist = t;
    statelist = m;

    for tsegment = 1:tsegments
        tsegment
        while t < tsegmentlength*tsegment
            %     if mod(round(t*100000)/100000,1)==0
            %         disp(num2str(t));
            %     end
            dt = random('Exponential',steptime(m),1,1);
            decider = random('Uniform',1E-10,1-1E-10,1,1);
            if decider <= Pfwd(m)
                m = m+1;
            else
                m = m-1;
            end
            t = t+dt;
            statelist = cat(1,statelist,m);
            tlist = cat(1,tlist,t);
            %     plot(tlist,statelist);
            %     pause(0.01);
        end

        if calculate_FRET == 1
            for k = 1:size(statelist)
                statelist(k) = FRETstates(statelist(k));
            end
        end

        if plot_average == 1
            % Perform downsampling to simulate camera acquisition
            disp('Performing downsampling...');
            bin = itime;
            k = 1;
            statelistbin = zeros(tsegmentlength/itime-1,1);
            currentbinstates = zeros(0,1);
            currentbintimes = zeros(0,1);
            currentbin_dt = zeros(0,1);
            while bin < tsegmentlength
                if mod(bin/itime,10)==0
                    disp(num2str(bin));
                end
                if tlist(k) <= bin
                    currentbinstates = cat(1,currentbinstates,statelist(k));
                    currentbintimes = cat(1,currentbintimes,tlist(k));
                    k = k+1;
                else
                    if mod(bin/itime,10)<0.5
                        disp('Done with...');
                        disp(num2str(bin));
                    end
                    currentbinstates = cat(1, currentbinstates, statelist(k));
                    currentbintimes = cat(1, currentbintimes, bin);
                    currentbintimes = currentbintimes - (bin - itime);
                    currentbintimes(length(currentbintimes)) = itime;
                    currentbin_dt(1,1) = currentbintimes(1);
                    for n = 2:length(currentbintimes)
                        currentbin_dt(n,1) = currentbintimes(n)-currentbintimes(n-1);
                    end
                    weights = (currentbin_dt./itime);
                    weightedlist = currentbinstates.*weights;
                    statelistbin(round(bin/itime)) = sum(weightedlist);
                    stateslistbin = zeros(tsegmentlength/itime);
                    bin = bin+itime;
                    currentbinstates = zeros(0,1);
                    currentbintimes = zeros(0,1);
                    currentbin_dt = zeros(0,1);
                end
            end
            % binlist = itime:itime:tfinal-itime;
            binlist = 1:1:length(stateslistbin);
            binlist = (binlist-1)*itime;
        end
        statelistbin_all = cat(1,statelistbin_all,statelistbin);
    end

    binlist_all = 1:1:length(statelistbin_all);
    binlist_all = (binlist_all-1)*itime;

    figure(1)
    plot(tlist,statelist,'k-','LineWidth',1);
    xlabel('Time (s)');
    ylabel('State');
    xlim([0 tfinal]);
    ylim([0 nstates-1]);
    % pause(0.01);
    set(gca,'FontSize', 24,'LineWidth',2);

    n_transitions = 0;
    for k = 2:length(statelist)-1
        if statelist(k) == b+3
            if statelist(k-1) ~= statelist(k+1)
                n_transitions = n_transitions+1;
            end
        end
    end
    transitions(run) = n_transitions;
    lifetime(run) = tfinal/n_transitions;

    if plot_average==1
        hold on
        plot(binlist,statelistbin,'r','LineWidth',2);
        % xlabel('Time (s)');
        % ylabel('Average State');
        hold off
    end

    traces = statelistbin';
    tracesname = strcat('W_a_',num2str(a),'_b_',num2str(b),'_kbindratio_',num2str(kbindratio),'_run_',num2str(run),'.dat');
    save(tracesname, 'traces');

end

% For autocorrelation:

% N = xcorr(matrix(:,2),matrix(:,2),'unbiased');
% len2 = (size(N,1)+1)/2;
% N1 = flipud(N(1:len2));
% N2 = N(len2:len2*2-1);
% Nav = (N1+N2)./2;
% Nav = Nav-min(Nav(1:round(len2/2)));
% Nav = Nav./max(Nav(1:round(len2/2)));
% figure(4)
% plot(0:1/kbm:(len2-1)*1/kbm,Nav,'k.-');
% ylabel('Norm. Autocorrelation');
% xlabel('Time (s)');
% plot(0:itime:itime*(length(Nav)-1),Nav,'k.-');
% tau = 0:itime:itime*(length(Nav)-1);
% cftool(tau(1:20),Nav(1:20));
