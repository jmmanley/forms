%% Test of forms code from
% Cashman, Thomas J., and Andrew W. Fitzgibbon. "What shape are dolphins? building 3d morphable models from 2d images." IEEE transactions on pattern analysis and machine intelligence 35.1 (2013): 232-244.
%
% Jason Manley, updated Mar 2018
% Question? jmanley@rockefeller.edu

addpath(genpath('/path/to/repo/forms/'))

%% Bananas

% read bananas project
bananas = read_project('projects/bananas.fpj');

% align images
bananas = align_project(bananas);

% set parameters
parameters = struct('M',2,'xi_0',0.5,'xi_def',0.25,'sigma_norm',4);
parameters = makeModelParameters(parameters);

% fit model
banana_model = forms(bananas, parameters);

% show image n before/after optimization
n=1;
figure;
subplot(1,2,1); title('Initial Conditions');
plot_modelfit(bananas, n);
subplot(1,2,2); title('Model Fit');
plot_modelfit(bananas, n, banana_model);

% plot the nth basis shape norm
n = 1;
figure;
plot_basisshapenorm(bananas, banana_model, n);

% plot all model fits
for i=1:length(bananas.images)
    plot_modelfit(bananas, i, banana_model);
    pause;
end


%% Dolphins

% read dolphins project
dolphins = read_project('projects/dolphins.fpj');

% align images
dolphins = align_project(dolphins);

% set parameters
parameters = struct('M',8,'xi_0',0.5,'xi_def',0.25,'sigma_norm',10);
parameters = makeModelParameters(parameters);

% fit model
dolphin_model = forms(dolphins, parameters);

% show image n before/after optimization
n=1;
figure;
subplot(1,2,1); title('Initial Conditions');
plot_modelfit(dolphins, n);
subplot(1,2,2); title('Model Fit');
plot_modelfit(dolphins, n, dolphin_model);

% plot the nth basis shape norm
n = 1;
figure;
plot_basisshapenorm(dolphins, dolphin_model, n);

% plot all model fits
figure;
for i=1:length(dolphins.images)
    plot_modelfit(dolphins, i, dolphin_model);
    pause;
end