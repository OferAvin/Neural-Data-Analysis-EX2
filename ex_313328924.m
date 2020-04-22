%% preparing and orgenizing date
clear all
close all
load SpikesX10U12D.mat
num_of_neurons = size(SpikesX10U12D,1);
num_of_degs = size(SpikesX10U12D,2);
num_of_rep = size(SpikesX10U12D,3);
bin_duration = 0.02;
times_bins_vec = (0:bin_duration:1.28);
chosen_neuron = 3;
num_of_plots_per_raws = 2;
num_of_bins = length(times_bins_vec);
mat = zeros(num_of_neurons,num_of_degs,num_of_rep,num_of_bins-1);
rate = zeros(num_of_degs,num_of_bins-1);
%% creating a 4dim array 
for i = 1:num_of_neurons
    for j = 1:num_of_degs
        for k = 1:num_of_rep
            mat(i,j,k,:) = histcounts(SpikesX10U12D(i,j,k).TimeList,times_bins_vec);
        end
    end
end
%% calculating rate for each degree for a chosen neuron
onesVec = ones(1,num_of_rep);
for i = 1:num_of_degs
    rep_bin_mat = squeeze(mat(chosen_neuron,i,:,:));
    sum_of_all_bins = onesVec*rep_bin_mat;
    rate(i,:) = sum_of_all_bins/(bin_duration*num_of_rep);
end
%% Ploting
figure();
hold on
for plotID = 1:num_of_degs
    subplot(num_of_plots_per_raws,num_of_degs/num_of_plots_per_raws,plotID);
    hold on
    if(plotID == 1 || plotID == 7)
        ylabel('rate[Hz]');
    end
    if(plotID >=7 && plotID <= 12)
        xlabel('time[sec]');
    end
    set(gca,'YLim',[0 30],'XLim',[0 1.28]);
    degree = num2str(plotID*30-30);
    title("\theta = "+ degree+" \circ");
    bar(times_bins_vec(1:end-1),rate(plotID,:),'histc');
    pbaspect([1 1 1])
end


  
