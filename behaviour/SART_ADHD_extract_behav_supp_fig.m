%% Init
clear all; 
close all;

if isempty(findstr(pwd,'thandrillon'))==0 
    path_LSCPtools='/Users/thandrillon/WorkGit/LSCPtools/';
    path_fieldtrip='/Users/thandrillon/WorkGit/projects/ext/fieldtrip/';
    behav_path = '/Users/thandrillon/Data/ADHD_MW/Behaviour/';
    save_path = '/Users/thandrillon/WorkGit/projects/inprogress/MWMB_ADHD/tables';
    path_RainCloudPlot='/Users/thandrillon/WorkGit/projects/ext/RainCloudPlots/';
else
    path_LSCPtools = '/Users/elaine/desktop/MATLAB_Functions/LSCPtools/';
    path_fieldtrip = '/Users/elaine/desktop/MATLAB_Functions/fieldtrip/';
    path_RainCloudPlot = '/Users/elaine/desktop/MATLAB_Functions/RainCloudPlots/';
    path_DataViz = '/Users/elaine/desktop/MATLAB_Functions/DataViz/';
    path_chi2test = '/Users/elaine/desktop/MATLAB_Functions/chi2test/';
    behav_path = '/Volumes/Seagate/MWMB_ADHD_SART/Behaviour/';
    preproc_path='/Volumes/Seagate/MWMB_ADHD_SART/preproc/';
    save_path = '/Users/elaine/desktop/Git_Studies/MWMB_ADHD/tables';
    
    %     mkdir(path_detectSW)
end
% adding relevant toolboxes to the path
addpath(genpath(path_LSCPtools))
addpath(genpath(path_RainCloudPlot));
addpath(genpath(path_DataViz));
addpath(path_fieldtrip)
addpath(path_chi2test)
ft_defaults;
% addpath(genpath(path_ExGauss))
% addpath(genpath(path_FMINSEARCHBND))

addpath(behav_path)

% select relevant files, here baseline blocks
files=dir([behav_path filesep 'wanderIM_behavres_*.mat']);

% NOTE: test_res columns: 
% nblock this_blockcond thiset ntrial this_seq_trial TargetID thisresp stimonset dur_face_no_rand thisresptime  this_nogo this_go

%NOTE: probe_res columns:
% this_probe this_probetime startProbe nblock this_blockcond ntrial probe_responses(:,1)' probe_responses(:,2)' probe_responses(:,3)' probe_responses(:,4)';

% Loading data 
all_behav_table = readtable([save_path filesep 'MWMB_ADHD_all_behav_byTrial.txt']);
all_probe_table = readtable([save_path filesep 'MWMB_ADHD_all_probe_behav.txt']);
all_block_table = readtable([save_path filesep 'MWMB_ADHD_all_block.txt']);

all_behav_table.Group=categorical(all_behav_table.Group);
all_behav_table.SubID=categorical(all_behav_table.SubID);

all_block_table.Group=categorical(all_block_table.Group);
all_block_table.SubID=categorical(all_block_table.SubID);

% Loading table with demo + questionnaire + byTrial data 
behav_demo_table = readtable([save_path filesep 'SART_ADHD_behav_demo_v1.txt']);
behav_demo_table.Group=categorical(behav_demo_table.Group);
behav_demo_table.SubID=categorical(behav_demo_table.SubID);
behav_demo_table.Sex=categorical(behav_demo_table.Sex);
behav_demo_table.ReportedSubtype=categorical(behav_demo_table.ReportedSubtype);
behav_demo_table.DIVASubtype=categorical(behav_demo_table.DIVASubtype);



%% Figures
Colors=[253,174,97;
    171,217,233;
    44,123,182]/256;

ctrs=unique(behav_demo_table.SubID(behav_demo_table.Group=='C' ));
Miss_CTR=[];
for nc=1:length(ctrs)
    Miss_CTR(nc) = nanmean(behav_demo_table.Misses(behav_demo_table.SubID == ctrs(nc)));
end

adhds_inn=unique(behav_demo_table.SubID(behav_demo_table.Group=='A' & behav_demo_table.ReportedSubtype =='Inattention'));
Miss_ADHD_inn=[];
for nc=1:length(adhds_inn)
    Miss_ADHD_inn(nc)=nanmean(behav_demo_table.Misses(behav_demo_table.SubID==adhds_inn(nc)));
end
Miss_ADHD_com = [];
adhds_com=unique(behav_demo_table.SubID(behav_demo_table.Group=='A' & behav_demo_table.ReportedSubtype =='Combined'));
for nc=1:length(adhds_com)
    Miss_ADHD_com(nc)=nanmean(behav_demo_table.Misses(behav_demo_table.SubID==adhds_com(nc)));
end

FA_CTR=[];
for nc=1:length(ctrs)
    FA_CTR(nc)=nanmean(all_behav_table.FA(all_behav_table.SubID==ctrs(nc)));
end
FA_ADHD_inn=[];
for nc=1:length(adhds_inn)
    FA_ADHD_inn(nc)=nanmean(all_behav_table.FA(all_behav_table.SubID==adhds_inn(nc)));
end
FA_ADHD_com=[];
for nc=1:length(adhds_com)
    FA_ADHD_com(nc)=nanmean(all_behav_table.FA(all_behav_table.SubID==adhds_com(nc)));
end

Hit_RT_CTR=[];
for nc=1:length(ctrs)
    Hit_RT_CTR(nc)=nanmean(all_behav_table.RT(all_behav_table.SubID==ctrs(nc)));
end
Hit_RT_ADHD_inn=[];
for nc=1:length(adhds_inn)
    Hit_RT_ADHD_inn(nc)=nanmean(all_behav_table.RT(all_behav_table.SubID==adhds_inn(nc)));
end
Hit_RT_ADHD_com=[];
for nc=1:length(adhds_com)
    Hit_RT_ADHD_com(nc)=nanmean(all_behav_table.RT(all_behav_table.SubID==adhds_com(nc)));
end

% EP - below I changed it to retrieve data from all_block_table instead of all_behav
stdRT_CTR=[];
for nc=1:length(ctrs)
    stdRT_CTR(nc)=nanmean(all_block_table.stdRT(all_block_table.SubID==ctrs(nc)));
end
stdRT_ADHD_inn=[];
for nc=1:length(adhds_inn)
    stdRT_ADHD_inn(nc)=nanmean(all_block_table.stdRT(all_block_table.SubID==adhds_inn(nc)));
end
stdRT_ADHD_com=[];
for nc=1:length(adhds_com)
    stdRT_ADHD_com(nc)=nanmean(all_block_table.stdRT(all_block_table.SubID==adhds_com(nc)));
end


%%
all_Miss=100*[Miss_CTR , Miss_ADHD_inn, Miss_ADHD_com]';
all_Group=[zeros(1,length(Miss_CTR)), ones(1,length(Miss_ADHD_inn)), 2 * ones(1,length(Miss_ADHD_com))]';

%figure; subplot(1,4,1)
figure('Position',[347   167   441   447]);
h = daviolinplot(all_Miss,'groups',all_Group,'outsymbol','k+','colors',Colors,...
    'boxcolors','same','scatter',1,'jitter',1,'xtlabels', {'Control','ADHD: Inattentive', 'ADHD: Combined'},...
    'boxwidth',1.5,'scattersize',70);
ylabel('% Misses');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
set(gca,'FontSize',10);
format_fig;
ylim([-0.5 15])

%%
all_FA=100*[FA_CTR , FA_ADHD_inn, FA_ADHD_com]';
all_Group=[zeros(1,length(FA_CTR)), ones(1,length(FA_ADHD_inn)), 2 * ones(1,length(FA_ADHD_com))]';

%subplot(1,4,2)
figure('Position',[347   167   441   447]);
h = daviolinplot(all_FA,'groups',all_Group,'outsymbol','k+','colors',Colors,...
    'boxcolors','same','scatter',1,'jitter',1,'xtlabels', {'Control','ADHD: Inattentive', 'ADHD: Combined'},...
    'boxwidth',1.5,'scattersize',70);
ylabel('% False Alarms');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
set(gca,'FontSize',10);
format_fig;
ylim([0 100])

%%
all_RT = [Hit_RT_CTR , Hit_RT_ADHD_inn, Hit_RT_ADHD_com]'; 
all_Group=[zeros(1,length(Hit_RT_CTR)), ones(1,length(Hit_RT_ADHD_inn)), 2 * ones(1,length(Hit_RT_ADHD_com))]';

%subplot(1,4,3)
figure('Position',[347   167   441   447]);
h = daviolinplot(all_RT,'groups',all_Group,'outsymbol','k+','colors',Colors,...
    'boxcolors','same','scatter',1,'jitter',1,'xtlabels', {'Control','ADHD: Inattentive', 'ADHD: Combined'},...
    'boxwidth',1.5,'scattersize',70);
ylabel('Reaction Time (s)');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
set(gca,'FontSize',10);
format_fig;
ylim([0.2 0.7])

%%
all_stdRT = [stdRT_CTR, stdRT_ADHD_inn stdRT_ADHD_com]'; 
all_Group=[zeros(1,length(stdRT_CTR)) , ones(1,length(stdRT_ADHD_inn)), 2 * ones(1,length(stdRT_ADHD_com))]';

%subplot(1,4,4)
figure('Position',[347   167   441   447]);
h = daviolinplot(all_stdRT,'groups',all_Group,'outsymbol','k+','colors',Colors,...
    'boxcolors','same','scatter',1,'jitter',1,'xtlabels', {'CTRL','ADHD: Inattentive', 'ADHD: Combined'},...
    'boxwidth',1.5,'scattersize',70);
ylabel('StdDev of Reaction Time (s)');
xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
set(gca,'FontSize',10);
format_fig;
ylim([0. 0.35])


%% Mind States
all_probe_table.SubID = string(all_probe_table.SubID); % Convert cell array to string
all_probe_table.Group = string(all_probe_table.Group);
behav_demo_table.SubID = string(behav_demo_table.SubID); % Convert character array to string
behav_demo_summary = varfun(@(x) x(1), behav_demo_table, ...
    'GroupingVariables', 'SubID', ...
    'InputVariables', {'ReportedSubtype'}); %Summary table because there were SubID repeats as it was byTrial
behav_demo_summary.Properties.VariableNames{'Fun_ReportedSubtype'} = 'ReportedSubtype';

probe_subtype_table = join(all_probe_table, behav_demo_summary(:, {'SubID', 'ReportedSubtype'}), 'Keys', 'SubID');
probe_subtype_table.SubID = categorical(probe_subtype_table.SubID); 
probe_subtype_table.Group = categorical(probe_subtype_table.Group); 


    numbers_of_interest = [1, 2, 3, 4];
    Control_state =[];
    ADHD_inn_state = [];
    ADHD_com_state = [];
    Control_state = zeros(length(ctrs), length(numbers_of_interest));
    ADHD_inn_state = zeros(length(adhds_inn), length(numbers_of_interest));
    ADHD_com_state = zeros(length(adhds_com), length(numbers_of_interest));

    clear state_values
    index = 1; % Initialize a separate index for Control_state
    for nc = 1:length(ctrs)
        if ctrs(nc) == 'C015'
            disp('Skipping C015')
            continue; % Skipping C015 'cause ppt didn't understand instructions
        end
        % Extract the State column for the current participant
        state_values = probe_subtype_table.State(probe_subtype_table.SubID == ctrs(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        Control_state(index, :) = histcounts(state_values, [numbers_of_interest, Inf]);
        index = index + 1; % Increment the Control_state index only when valid data is added
    end

    % Calculate the total occurrences of each number across all participants
    Ctr_state_total_occurrences = sum(Control_state);
    % Calculate the percentage distribution
    Ctr_state_percentage_distribution = Ctr_state_total_occurrences / sum(Ctr_state_total_occurrences) * 100;


    clear state_values
    for nc = 1:length(adhds_inn)
        % Extract the State column for the current participant
        state_values = probe_subtype_table.State(probe_subtype_table.SubID == adhds_inn(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_inn_state(nc, :) = histcounts(state_values, [numbers_of_interest, Inf]);
    end

    ADHD_inn_state_total_occurrences = sum(ADHD_inn_state);
    ADHD_inn_state_percentage_distribution = ADHD_inn_state_total_occurrences / sum(ADHD_inn_state_total_occurrences) * 100;

    clear state_values
    for nc = 1:length(adhds_com)
        % Extract the State column for the current participant
        state_values = probe_subtype_table.State(probe_subtype_table.SubID == adhds_com(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_com_state(nc, :) = histcounts(state_values, [numbers_of_interest, Inf]);
    end

    ADHD_com_state_total_occurrences = sum(ADHD_com_state);
    ADHD_com_state_percentage_distribution = ADHD_com_state_total_occurrences / sum(ADHD_com_state_total_occurrences) * 100;



    % Plot the distribution - bar graph
    labels = {'On Task', 'Mind Wandering','Mind Blanking', 'Dont Remember'};
    numGroups = 3;  % Control, ADHD_inn, ADHD_com
    groupOffsets = linspace(-0.2, 0.2, numGroups);  % Calculate offsets for each group
    barWidth = 0.29;

    figure('Position',[1955         361         758         413]);
    barPositions = (1:numel(labels)) - 0.3;  % Adjust position for Control group
    h1 = bar(barPositions, Ctr_state_percentage_distribution, 'FaceColor', Colors(1,:), 'BarWidth', barWidth);

    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(Ctr_state_percentage_distribution, 1)
            text(barPositions(i), Ctr_state_percentage_distribution(j, i), sprintf('%.1f%%', Ctr_state_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25,'Color',Colors(1,:),'FontWeight','bold');
        end
    end
    format_fig

    hold on;
    barPositions = (1:numel(labels)) + 0.0;  % Adjust position for ADHD Inattention group
    h2 = bar(barPositions, ADHD_inn_state_percentage_distribution, 'FaceColor', Colors(2,:), 'BarWidth', barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(ADHD_inn_state_percentage_distribution, 1)
            text(barPositions(i), ADHD_inn_state_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_inn_state_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25,'Color',Colors(2,:),'FontWeight','bold');
        end
    end
    format_fig
    set(gca,'FontSize',30,'FontWeight','bold')
    xlim([0.5 4.5])

    barPositions = (1:numel(labels)) + 0.3;  % Adjust position for ADHD Combined group
    h3 = bar(barPositions, ADHD_com_state_percentage_distribution, 'FaceColor', Colors(3,:), 'BarWidth', barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(ADHD_com_state_percentage_distribution, 1)
            text(barPositions(i), ADHD_com_state_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_com_state_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25,'Color',Colors(3,:),'FontWeight','bold');
        end
    end

    % Set x-axis labels and ticks
    xticks(1:numel(labels));
    xticklabels(labels);
    xtickangle(45);
    xlim([0.5 4.5]);
    ylabel('% of Mind State');
    ylim([0 70]);
    set(gca, 'FontSize', 30, 'FontWeight', 'bold');
    legend([h1, h2, h3], {'Control', 'ADHD: Inattentive', 'ADHD: Combined'}, 'Location', 'northeast', 'FontSize', 18);
    format_fig;


    % Plot the distribution - raincloud plots
    labels = {'ON', 'MW','MB', '??'};
    all_States=[Control_state ; ADHD_inn_state ; ADHD_com_state]./40*100;
    all_Group=[zeros(size(Control_state,1),1) ; ones(size(ADHD_inn_state,1),1) ; 2 * ones(size(ADHD_com_state,1),1)]';

    figure('Position',[347         335        1028         279]);
    h = daviolinplot(all_States,'groups',all_Group,'outsymbol','k+','colors',Colors,...
        'boxcolors','same','scatter',1,'jitter',1,'xtlabels', labels,...
        'boxwidth',1,'scattersize',30);
    ylabel('% Responses');
    xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
    set(gca,'FontSize',10);
    format_fig;
    ylim([0 100])


 %% Mind Wandering by Distraction
    numbers_of_interest = [1, 2, 3];
    ADHD_inn_mw = [];
    ADHD_com_mw = [];
    Control_mw =[];
    Control_mw = zeros(length(ctrs), length(numbers_of_interest));
    ADHD_inn_mw = zeros(length(adhds_inn), length(numbers_of_interest));
    ADHD_com_mw= zeros(length(adhds_com), length(numbers_of_interest));

    clear mw_values
    index = 1; % Initialize a separate index for Control_vig
    for nc = 1:length(ctrs)
        if ctrs(nc) == 'C015'
            disp('Skipping C015')
            continue; % Skipping C015 'cause ppt didn't understand instructions
        end
        % Extract the State column for the current participant
        mw_values = probe_subtype_table.Distraction(probe_subtype_table.SubID == ctrs(nc));
        % Count occurrences of each number (1, 2, 3) for the current participant
        Control_mw(index, :) = histcounts(mw_values, [numbers_of_interest, Inf]);
        index = index + 1; % Increment the Control_state index only when valid data is added
    end

    % Calculate the total occurrences of each number across all participants
    Ctr_mw_total_occurrences = sum(Control_mw);
    % Calculate the percentage distribution
    Ctr_mw_percentage_distribution = Ctr_mw_total_occurrences / sum(Ctr_mw_total_occurrences) * 100;

    clear mw_values
    for nc = 1:length(adhds_inn)
        % Extract the State column for the current participant
        mw_values = probe_subtype_table.Distraction(probe_subtype_table.SubID == adhds_inn(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_inn_mw(nc, :) = histcounts(mw_values, [numbers_of_interest, Inf]);
    end

    ADHD_inn_mw_total_occurrences = sum(ADHD_inn_mw);
    ADHD_inn_mw_percentage_distribution = ADHD_inn_mw_total_occurrences / sum(ADHD_inn_mw_total_occurrences) * 100;

    clear mw_values
    for nc = 1:length(adhds_com)
        % Extract the State column for the current participant
        mw_values = probe_subtype_table.Distraction(probe_subtype_table.SubID == adhds_com(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_com_mw(nc, :) = histcounts(mw_values, [numbers_of_interest, Inf]);
    end

    ADHD_com_mw_total_occurrences = sum(ADHD_com_mw);
    ADHD_com_mw_percentage_distribution = ADHD_com_mw_total_occurrences / sum(ADHD_com_mw_total_occurrences) * 100;



    % Plot the distribution - bar graphs
    labels = {'Smth in Room', 'Personal','About task'};
    numGroups = 3;  % Control, ADHD_inn, ADHD_com
    groupOffsets = linspace(-0.2, 0.2, numGroups);  % Calculate offsets for each group
    barWidth = 0.29;
    figure%('Position',[100         100         800         600]);

    barPositions = (1:numel(labels)) - 0.3; 
    h1 = bar(barPositions, Ctr_mw_percentage_distribution, 'FaceColor',Colors(1,:),'BarWidth',barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(Ctr_mw_percentage_distribution, 1)
            text(barPositions(i), Ctr_mw_percentage_distribution(j, i), sprintf('%.1f%%', Ctr_mw_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom','FontSize', 25,'Color',Colors(1,:),'FontWeight','bold');
        end
    end
    hold on;

    barPositions = (1:numel(labels)) + 0.0;
    h2 = bar(barPositions, ADHD_inn_mw_percentage_distribution, 'FaceColor',Colors(2,:),'BarWidth',barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(ADHD_inn_mw_percentage_distribution, 1)
            text(barPositions(i), ADHD_inn_mw_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_inn_mw_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25,'Color',Colors(2,:),'FontWeight','bold');
        end
    end

    barPositions = (1:numel(labels)) + 0.3;
    h3 = bar(barPositions, ADHD_com_mw_percentage_distribution, 'FaceColor',Colors(3,:),'BarWidth',barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(ADHD_com_mw_percentage_distribution, 1)
            text(barPositions(i), ADHD_com_mw_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_com_mw_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25,'Color',Colors(3,:),'FontWeight','bold');
        end
    end
    % Set x-axis labels and ticks
    xticks(1:numel(labels));
    xticklabels(labels);
    xtickangle(45);
    xlim([0.5 numel(labels) + 0.5]);
    ylabel('% of Types of MW');
    ylim([0 80]);
    set(gca, 'FontSize', 30, 'FontWeight', 'bold');
    legend([h1, h2, h3], {'Control', 'ADHD: Inattentive', 'ADHD: Combined'}, 'Location', 'northeast', 'FontSize', 18);
    format_fig;


    %% Distraction - Pair Wise Chi squared tests
distraction_levels = 1:3;

% Initialize result arrays
p_values = NaN(length(distraction_levels), 3); % 3 comparisons (Inattention vs. Control, Combined vs. Control, Inattention vs. Combined)
Q_values = NaN(length(distraction_levels), 3);

% Exclude rows with NaN in Distraction
valid_rows = ~isnan(probe_subtype_table.Distraction);

% Filter the table to only valid rows
filtered_table = probe_subtype_table(valid_rows, :);

% Loop through each distraction level
for distraction_level = distraction_levels
    % Get counts for each subtype
    adhd_inn_distraction = sum(filtered_table.Distraction == distraction_level & filtered_table.ReportedSubtype == 'Inattention');
    adhd_com_distraction = sum(filtered_table.Distraction == distraction_level & filtered_table.ReportedSubtype == 'Combined');
    control_distraction = sum(filtered_table.Distraction == distraction_level & filtered_table.ReportedSubtype == 'Control');

    % Get counts for each subtype not selecting the current distraction level
    adhd_inn_not_distraction = sum(filtered_table.Distraction ~= distraction_level & filtered_table.ReportedSubtype == 'Inattention');
    adhd_com_not_distraction = sum(filtered_table.Distraction ~= distraction_level & filtered_table.ReportedSubtype == 'Combined');
    control_not_distraction = sum(filtered_table.Distraction ~= distraction_level & filtered_table.ReportedSubtype == 'Control');

    % Pairwise comparisons
    % Inattention vs. Control
    count_matrix_inn_control = [adhd_inn_distraction, adhd_inn_not_distraction;
                                control_distraction, control_not_distraction];
    [p1, Q1] = chi2test(count_matrix_inn_control);

    % Combined vs. Control
    count_matrix_com_control = [adhd_com_distraction, adhd_com_not_distraction;
                                control_distraction, control_not_distraction];
    [p2, Q2] = chi2test(count_matrix_com_control);

    % Inattention vs. Combined
    count_matrix_inn_com = [adhd_inn_distraction, adhd_inn_not_distraction;
                            adhd_com_distraction, adhd_com_not_distraction];
    [p3, Q3] = chi2test(count_matrix_inn_com);

    % Store results
    p_values(distraction_level, :) = [p1, p2, p3];
    Q_values(distraction_level, :) = [Q1, Q2, Q3];

    % Display results for this distraction level
    disp(['Distraction Level ', num2str(distraction_level)]);
    disp('Inattention vs. Control:');
    disp(count_matrix_inn_control);
    disp(['p-value: ', num2str(p1), ', Q-statistic: ', num2str(Q1)]);

    disp('Combined vs. Control:');
    disp(count_matrix_com_control);
    disp(['p-value: ', num2str(p2), ', Q-statistic: ', num2str(Q2)]);

    disp('Inattention vs. Combined:');
    disp(count_matrix_inn_com);
    disp(['p-value: ', num2str(p3), ', Q-statistic: ', num2str(Q3)]);
    disp('------------------------------');
end

% Display summary of results
disp('Summary of Chi-square test results for all distraction levels:');
summary_table = table(distraction_levels', p_values, Q_values, ...
                      'VariableNames', {'DistractionLevel', 'pValues', 'QStatistics'});
disp(summary_table);

% Bonferroni correction
num_tests = 3; % Number of pairwise comparisons per distraction level
corrected_alpha = 0.05 / num_tests;

% Check for significance after Bonferroni correction
significant_results = p_values < corrected_alpha;

% Display corrected results
disp('Summary with Bonferroni correction:');
for distraction_level = 1:length(distraction_levels)
    disp(['Distraction Level ', num2str(distraction_level)]);
    disp(['Corrected alpha: ', num2str(corrected_alpha)]);
    disp(['Significant Comparisons (Inattention vs. Control, Combined vs. Control, Inattention vs. Combined): ', ...
        num2str(significant_results(distraction_level, :))]);
end


%% Intentionality  (Note: 1 = Entirely intentional to 4 = Entirely unintentional)
  numbers_of_interest = [1, 2, 3, 4];
    ADHD_inn_int = [];
    ADHD_com_int = [];
    Control_int =[];
    Control_int = zeros(length(ctrs), length(numbers_of_interest));
    ADHD_inn_int = zeros(length(adhds_inn), length(numbers_of_interest));
    ADHD_com_int= zeros(length(adhds_com), length(numbers_of_interest));

    clear int_values
    index = 1; % Initialize a separate index for Control_vig
    for nc = 1:length(ctrs)
        if ctrs(nc) == 'C015'
            disp('Skipping C015')
            continue; % Skipping C015 'cause ppt didn't understand instructions
        end
        % Extract the State column for the current participant
        int_values = probe_subtype_table.Intention(probe_subtype_table.SubID == ctrs(nc) & probe_subtype_table.StateC == 'MW'); %EP - changed this so that it only counts intention values when they report MW
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        Control_int(index, :) = histcounts(int_values, [numbers_of_interest, Inf]);
        index = index + 1; % Increment the Control_int index only when valid data is added
    end

    % Calculate the total occurrences of each number across all participants
    Ctr_int_total_occurrences = sum(Control_int);
    % Calculate the percentage distribution
    Ctr_int_percentage_distribution = Ctr_int_total_occurrences / sum(Ctr_int_total_occurrences) * 100;

    clear int_values
    for nc = 1:length(adhds_inn)
        % Extract the State column for the current participant
        int_values = probe_subtype_table.Intention(probe_subtype_table.SubID == adhds_inn(nc) & probe_subtype_table.StateC == 'MW'); %EP - changed this so that it only counts intention values when they report MW
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_inn_int(nc, :) = histcounts(int_values, [numbers_of_interest, Inf]);
    end

    ADHD_inn_int_total_occurrences = sum(ADHD_inn_int);
    ADHD_inn_int_percentage_distribution = ADHD_inn_int_total_occurrences / sum(ADHD_inn_int_total_occurrences) * 100;

    clear int_values
    for nc = 1:length(adhds_com)
        % Extract the State column for the current participant
        int_values = probe_subtype_table.Distraction(probe_subtype_table.SubID == adhds_com(nc) & probe_subtype_table.StateC == 'MW'); %EP - changed this so that it only counts intention values when they report MW
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_com_int(nc, :) = histcounts(int_values, [numbers_of_interest, Inf]);
    end

    ADHD_com_int_total_occurrences = sum(ADHD_com_int);
    ADHD_com_int_percentage_distribution = ADHD_com_int_total_occurrences / sum(ADHD_com_int_total_occurrences) * 100;



    % Plot the distribution - bar graphs
    labels = {'Entirely Int.', 'Somewhat Int.','Somewhat Unint.', 'Entirely Unint.'};
    numGroups = 3;  % Control, ADHD_inn, ADHD_com
    groupOffsets = linspace(-0.2, 0.2, numGroups);  % Calculate offsets for each group
    barWidth = 0.29;
    figure%('Position',[100         100         800         600]);

    barPositions = (1:numel(labels)) - 0.3; 
    h1 = bar(barPositions, Ctr_int_percentage_distribution, 'FaceColor',Colors(1,:),'BarWidth',barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(Ctr_int_percentage_distribution, 1)
            text(barPositions(i), Ctr_int_percentage_distribution(j, i), sprintf('%.1f%%', Ctr_int_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom','FontSize', 25,'Color',Colors(1,:),'FontWeight','bold');
        end
    end
    hold on;

    barPositions = (1:numel(labels)) + 0.0;
    h2 = bar(barPositions, ADHD_inn_int_percentage_distribution, 'FaceColor',Colors(2,:),'BarWidth',barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(ADHD_inn_int_percentage_distribution, 1)
            text(barPositions(i), ADHD_inn_int_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_inn_int_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25,'Color',Colors(2,:),'FontWeight','bold');
        end
    end

    barPositions = (1:numel(labels)) + 0.3;
    h3 = bar(barPositions, ADHD_com_int_percentage_distribution, 'FaceColor',Colors(3,:),'BarWidth',barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(ADHD_com_int_percentage_distribution, 1)
            text(barPositions(i), ADHD_com_int_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_com_int_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25,'Color',Colors(3,:),'FontWeight','bold');
        end
    end
    % Set x-axis labels and ticks
    xticks(1:numel(labels));
    xticklabels(labels);
    xtickangle(45);
    xlim([0.5 numel(labels) + 0.5]);
    ylabel('% of Intention Responses for MW');
    ylim([0 80]);
    set(gca, 'FontSize', 30, 'FontWeight', 'bold');
    legend([h1, h2, h3], {'Control', 'ADHD: Inattentive', 'ADHD: Combined'}, 'Location', 'northeast', 'FontSize', 18);
    format_fig;


    probe_subtype_table.StateC=categorical(nan(size(probe_subtype_table,1),1));
    probe_subtype_table.StateC(probe_subtype_table.State==1)='ON';
    probe_subtype_table.StateC(probe_subtype_table.State==2)='MW';
    probe_subtype_table.StateC(probe_subtype_table.State==3)='MB';
    probe_subtype_table.StateC(probe_subtype_table.State==4)='DK';

    probe_subtype_table.ReportedSubtype((probe_subtype_table.Group== 'C')) = {'Control'}; 
    
    figure;
    probe_subtype_table.ReportedSubtype = reordercats(probe_subtype_table.ReportedSubtype, {'Control', 'Inattention', 'Combined'});
    myReportedSubtype=categories(probe_subtype_table.ReportedSubtype);
    myStates=unique(probe_subtype_table.StateC);
    hb=[];
    Int_Paired_test=[];
    for nSta=1:4
        for_paired_states={};
        for nG=1:numel(myReportedSubtype)
            % Extract intention values for the current group and state
            state_values = probe_subtype_table.Intention(probe_subtype_table.StateC == myStates(nSta) & probe_subtype_table.ReportedSubtype == myReportedSubtype(nG));
            % Calculate percentage for current state and group
            state_percentages(nG, nSta) = numel(state_values) / numel(probe_subtype_table.Intention(probe_subtype_table.ReportedSubtype == myReportedSubtype(nG))) * 100;

            temp_bar=grpstats(probe_subtype_table.Intention(probe_subtype_table.StateC==myStates(nSta) & probe_subtype_table.ReportedSubtype==myReportedSubtype(nG)),...
                probe_subtype_table.SubID(probe_subtype_table.StateC==myStates(nSta) & probe_subtype_table.ReportedSubtype==myReportedSubtype(nG)));
            hb(nG)=simpleBarPlot(nSta+(1.5*nG-3)*0.2,temp_bar,Colors(nG,:),0.28,'none',[],2); %none was 'k' before to show SE (or was it SD?) bars

            % Add percentage text on top of each bar
            text(nSta + (1.5 * nG - 3) * 0.2, nanmean(temp_bar) + 0.01, sprintf('%.1f%%', state_percentages(nG, nSta)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20, 'FontWeight', 'bold', 'Color', Colors(nG,:));


            for_paired_states{nG}=temp_bar;
        end
        [p, tbl, stats] = anova1([for_paired_states{1}; for_paired_states{2}; for_paired_states{3}], ...
            [ones(size(for_paired_states{1})); 2*ones(size(for_paired_states{2})); 3*ones(size(for_paired_states{3}))],'off');

        %[h,p,ci,stats] = ttest2(for_paired_states{1},for_paired_states{2}); % replaced t test as now we're comparing 3 groups
%         Int_Paired_test(nSta,1)=p;
%         Int_Paired_test(nSta,2)=stats.tstat;
%         Int_Paired_test(nSta,3)=stats.df;
    end
    xlabels = {'On Task', 'Mind Wandering','Mind Blanking', 'Dont Remember'};
    ylabels = {'Entirely Int.', 'Somewhat Int.','Somewhat Unint.', 'Entirely Unint.'};
    set(gca,'Xtick',1:4,'XTickLabel',xlabels);
    set(gca,'Ytick',1:4,'YTickLabel',ylabels);
    ytickangle(45);
    xtickangle(45);
    legend(hb,myReportedSubtype)
    format_fig;
    ylabel('Intentionality')
    %% Vigilance
    numbers_of_interest = [1, 2, 3, 4];
    Control_vig =[];
    ADHD_inn_vig = [];
    ADHD_com_vig = [];
    Control_vig = zeros(length(ctrs), length(numbers_of_interest));
    ADHD_inn_vig = zeros(length(adhds_inn), length(numbers_of_interest));
    ADHD_com_vig = zeros(length(adhds_com), length(numbers_of_interest));

    clear vig_values
    index = 1; % Initialize a separate index for Control_vig
    for nc = 1:length(ctrs)
        if ctrs(nc) == 'C015'
            disp('Skipping C015')
            continue; % Skipping C015 'cause ppt didn't understand instructions
        end
        % Extract the State column for the current participant
        vig_values = probe_subtype_table.Vigilance(probe_subtype_table.SubID == ctrs(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        Control_vig(index, :) = histcounts(vig_values, [numbers_of_interest, Inf]);
        index = index + 1; % Increment the Control_state index only when valid data is added
    end

    % Calculate the total occurrences of each number across all participants
    Ctr_vig_total_occurrences = sum(Control_vig);
    % Calculate the percentage distribution
    Ctr_vig_percentage_distribution = Ctr_vig_total_occurrences / sum(Ctr_vig_total_occurrences) * 100;

    clear vig_values
    for nc = 1:length(adhds_inn)
        % Extract the State column for the current participant
        vig_values = probe_subtype_table.Vigilance(probe_subtype_table.SubID == adhds_inn(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_inn_vig(nc, :) = histcounts(vig_values, [numbers_of_interest, Inf]);
    end

    ADHD_inn_vig_total_occurrences = sum(ADHD_inn_vig);
    ADHD_inn_vig_percentage_distribution = ADHD_inn_vig_total_occurrences / sum(ADHD_inn_vig_total_occurrences) * 100;

    clear vig_values
    for nc = 1:length(adhds_com)
        % Extract the State column for the current participant
        vig_values = probe_subtype_table.Vigilance(probe_subtype_table.SubID == adhds_com(nc));
        % Count occurrences of each number (1, 2, 3, 4) for the current participant
        ADHD_com_vig(nc, :) = histcounts(vig_values, [numbers_of_interest, Inf]);
    end

    ADHD_com_vig_total_occurrences = sum(ADHD_com_vig);
    ADHD_com_vig_percentage_distribution = ADHD_com_vig_total_occurrences / sum(ADHD_com_vig_total_occurrences) * 100;

    all_Vigs=[Control_vig ; ADHD_inn_vig ; ADHD_com_vig]./40*100;


    % Plot the distribution - bar graphs
    labels = {'Alert ++', 'Alert +','Sleepy +', 'Sleepy ++'};
    numGroups = 3;  % Control, ADHD_inn, ADHD_com
    groupOffsets = linspace(-0.2, 0.2, numGroups);  % Calculate offsets for each group
    barWidth = 0.29;

    figure('Position',[1955         361         758         413]);
    barPositions = (1:numel(labels)) - 0.3;  % Adjust position for Control group
    h1 = bar(barPositions, Ctr_vig_percentage_distribution, 'FaceColor', Colors(1,:), 'BarWidth', barWidth);

    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(Ctr_vig_percentage_distribution, 1)
            text(barPositions(i), Ctr_vig_percentage_distribution(j, i), sprintf('%.1f%%', Ctr_vig_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25,'Color',Colors(1,:),'FontWeight','bold');
        end
    end
    format_fig

    hold on;
    barPositions = (1:numel(labels)) + 0.0;  % Adjust position for ADHD Inattention group
    h2 = bar(barPositions, ADHD_inn_vig_percentage_distribution, 'FaceColor', Colors(2,:), 'BarWidth', barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(ADHD_inn_vig_percentage_distribution, 1)
            text(barPositions(i), ADHD_inn_vig_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_inn_vig_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25,'Color',Colors(2,:),'FontWeight','bold');
        end
    end
    format_fig
    set(gca,'FontSize',30,'FontWeight','bold')
    xlim([0.5 4.5])

    barPositions = (1:numel(labels)) + 0.3;  % Adjust position for ADHD Combined group
    h3 = bar(barPositions, ADHD_com_vig_percentage_distribution, 'FaceColor', Colors(3,:), 'BarWidth', barWidth);
    % Add percentage values on top of each bar
    for i = 1:numel(labels)
        for j = 1:size(ADHD_com_vig_percentage_distribution, 1)
            text(barPositions(i), ADHD_com_vig_percentage_distribution(j, i), sprintf('%.1f%%', ADHD_com_vig_percentage_distribution(j, i)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 25,'Color',Colors(3,:),'FontWeight','bold');
        end
    end

    % Set x-axis labels and ticks
    xticks(1:numel(labels));
    xticklabels(labels);
    xtickangle(45);
    xlim([0.5 4.5]);
    ylabel('% of Responses');
    ylim([0 70]);
    set(gca, 'FontSize', 30, 'FontWeight', 'bold');
    legend([h1, h2, h3], {'Control', 'ADHD: Inattentive', 'ADHD: Combined'}, 'Location', 'northeast', 'FontSize', 18);
    format_fig;


    % Plot the distribution - raincloud plots
    labels = {'Alert ++', 'Alert +','Sleepy +', 'Sleepy ++'};

    figure('Position',[347         335        1028         279]);
    h = daviolinplot(all_Vigs,'groups',all_Group,'outsymbol','k+','colors',Colors,...
        'boxcolors','same','scatter',1,'jitter',1,'xtlabels', labels,...
        'boxwidth',1,'scattersize',30);
    ylabel('% Responses');
    xl = xlim; xlim([xl(1)-0.1, xl(2)+0.2]); % make more space for the legend
    set(gca,'FontSize',10);
    format_fig;
    ylim([0 100])

%%  Prepping for stats
% Adding 'control' as a subtype for controls 
probe_subtype_table1 = probe_subtype_table;
probe_subtype_table1.ReportedSubtype((probe_subtype_table1.Group== 'C')) = {'Control'}; 
% Set 'Control' as the reference group
probe_subtype_table1.ReportedSubtype = reordercats(probe_subtype_table1.ReportedSubtype, {'Control', 'Combined', 'Inattention'});
probe_subtype_table1.ReportedSubtype = categorical(probe_subtype_table1.ReportedSubtype);
%probe_subtype_table1.State = categorical(probe_subtype_table1.State, [1, 2, 3, 4], {'On Task', 'Mind Wandering', 'Mind Blanking', 'Dont Remember'});;

behav_demo_summary.SubID = categorical(behav_demo_summary.SubID);
block_subtype_table = join(all_block_table, behav_demo_summary(:, {'SubID', 'ReportedSubtype'}), 'Keys', 'SubID');
block_subtype_table.SubID = categorical(block_subtype_table.SubID); 
block_subtype_table.Group = categorical(block_subtype_table.Group); 
block_subtype_table.ReportedSubtype((block_subtype_table.Group== 'C')) = {'Control'}; 
block_subtype_table.ReportedSubtype = reordercats(block_subtype_table.ReportedSubtype, {'Control', 'Combined', 'Inattention'});
block_subtype_table.ReportedSubtype = categorical(block_subtype_table.ReportedSubtype);

behav_demo_table1=behav_demo_table;
behav_demo_table1.ReportedSubtype((behav_demo_table1.Group== 'C')) = {'Control'}; 
behav_demo_table1.ReportedSubtype = reordercats(behav_demo_table1.ReportedSubtype, {'Control', 'Combined', 'Inattention'});
%% Behaviour (trial level data)
% FAs/Commission Errors
% mdlFA0  = fitglme(behav_demo_table1,'FA~1+BlockN+ReportedSubtype+(1|SubID)');
% mdlFA1  = fitglme(behav_demo_table1,'FA~1+BlockN*ReportedSubtype+(1|SubID)'); 
mdlFA2  = fitglme(behav_demo_table1,'FA~1+BlockN*ReportedSubtype+(BlockN|SubID)'); % Winning AIC and BIC model - Group: p = .02, BlockN p <.001
%  %%% Extract fit statistics for each model
% AIC_values = [mdlFA0.ModelCriterion.AIC, mdlFA1.ModelCriterion.AIC, mdlFA2.ModelCriterion.AIC];
% BIC_values = [mdlFA0.ModelCriterion.BIC, mdlFA1.ModelCriterion.BIC, mdlFA2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlFA2)

% Misses/Omission Errors
% mdlMiss0  = fitglme(behav_demo_table1,'Misses~1+BlockN+ReportedSubtype+(1|SubID)'); 
% mdlMiss1  = fitglme(behav_demo_table1,'Misses~1+BlockN*ReportedSubtype+(1|SubID)'); 
mdlMiss2  = fitglme(behav_demo_table1,'Misses~1+BlockN*ReportedSubtype+(BlockN|SubID)');
%  %%% Extract fit statistics for each model
% AIC_values = [mdlMiss0.ModelCriterion.AIC, mdlMiss1.ModelCriterion.AIC, mdlMiss2.ModelCriterion.AIC];
% BIC_values = [mdlMiss0.ModelCriterion.BIC, mdlMiss1.ModelCriterion.BIC, mdlMiss2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);
% % Another way to compare models: 
% compare(mdlMiss0, mdlMiss1)
% compare(mdlMiss1,mdlMiss2)
% compare(mdlMiss0,mdlMiss2)
anova(mdlMiss2)

%RT
mdlRT0  = fitlme(behav_demo_table1,'RT~1+BlockN+ReportedSubtype+(1|SubID)'); 
mdlRT1  = fitlme(behav_demo_table1,'RT~1+BlockN*ReportedSubtype+(1|SubID)'); 
mdlRT2  = fitlme(behav_demo_table1,'RT~1+BlockN*ReportedSubtype+(BlockN|SubID)'); % Winning model 
%  %% Extract fit statistics for each model
% AIC_values = [mdlRT0.ModelCriterion.AIC, mdlRT1.ModelCriterion.AIC, mdlRT2.ModelCriterion.AIC];
% BIC_values = [mdlRT0.ModelCriterion.BIC, mdlRT1.ModelCriterion.BIC, mdlRT2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlRT2)


% Standard Deviation of Reaction Times
% NOTE: EP changed to block_subtype_table from behav_demo_table1 as the latter table just has the stdRT duplicated over trials
% mdlstdRT0  = fitlme(block_subtype_table,'stdRT~1+BlockN+ReportedSubtype+(1|SubID)'); 
% mdlstdRT1  = fitlme(block_subtype_table,'stdRT~1+BlockN*ReportedSubtype+(1|SubID)'); 
mdlstdRT2  = fitlme(block_subtype_table,'stdRT~1+BlockN*ReportedSubtype+(BlockN|SubID)'); % Winning model
%%% Extract fit statistics for each model
% AIC_values = [mdlstdRT0.ModelCriterion.AIC, mdlstdRT1.ModelCriterion.AIC, mdlstdRT2.ModelCriterion.AIC];
% BIC_values = [mdlstdRT0.ModelCriterion.BIC, mdlstdRT1.ModelCriterion.BIC, mdlstdRT2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlstdRT2)
%EP - do a post hoc 


%% Mind states - (1 = On, 2 = MW, 3 = MB, 4 = DK)
% mdlstate0 = fitlme(probe_subtype_table1,'State~1+Block+ReportedSubtype+(1|SubID)');
% mdlstate1 = fitlme(probe_subtype_table1,'State~1+Block*ReportedSubtype+(1|SubID)');
mdlstate2 = fitlme(probe_subtype_table1,'State~1+Block*ReportedSubtype+(Block|SubID)'); % Winning AIC & BIC model; Group: p <.001, Block: p =.032
% %% Extract fit statistics for each model
% AIC_values = [mdlstate0.ModelCriterion.AIC, mdlstate1.ModelCriterion.AIC, mdlstate2.ModelCriterion.AIC];
% BIC_values = [mdlstate0.ModelCriterion.BIC, mdlstate1.ModelCriterion.BIC, mdlstate2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlstate2)

% ON_table = probe_subtype_table1(probe_subtype_table1.State == 1, :); % NOTE: I commented this out as we're underpowered to do it at the probe level
% mdlON = fitlme(ON_table, 'State~1+Block+ReportedSubtype+(Block|SubID)');

mdlON = fitlme(block_subtype_table,'ON~1+BlockN+ReportedSubtype+(1|SubID)');
anova(mdlON)

mdlMW=fitlme(block_subtype_table,'MW~1+BlockN+ReportedSubtype+(1|SubID)');
anova(mdlMW)

mdlMB=fitlme(block_subtype_table,'MB~1+BlockN+ReportedSubtype+(1|SubID)');
anova(mdlMB)

mdlDK=fitlme(block_subtype_table,'DK~1+BlockN+ReportedSubtype+(1|SubID)');
anova(mdlDK)

%% Vigilance (1:Alert++ to 4:Sleepy ++)
% mdlVig0=fitlme(probe_subtype_table1,'Vigilance~1+Block+ReportedSubtype+(1|SubID)'); 
% mdlVig1=fitlme(probe_subtype_table1,'Vigilance~1+Block*ReportedSubtype+(1|SubID)');
mdlVig2=fitlme(probe_subtype_table1,'Vigilance~1+Block*ReportedSubtype+(Block|SubID)');% Winning AIC & BIC Model
% %% Extract fit statistics for each model
% AIC_values = [mdlVig0.ModelCriterion.AIC, mdlVig1.ModelCriterion.AIC, mdlVig2.ModelCriterion.AIC];
% BIC_values = [mdlVig0.ModelCriterion.BIC, mdlVig1.ModelCriterion.BIC, mdlVig2.ModelCriterion.BIC];
% % Display results in a table
% ModelNames = {'Model 0', 'Model 1', 'Model 2'};
% fit_table = table(ModelNames', AIC_values', BIC_values', 'VariableNames', {'Model', 'AIC', 'BIC'});
% disp(fit_table);

anova(mdlVig2)

results_subtype = multcompare(anova(mdlVig2), 'ReportedSubtype', 'DFMethod', 'Satterthwaite');
disp(results_subtype);

