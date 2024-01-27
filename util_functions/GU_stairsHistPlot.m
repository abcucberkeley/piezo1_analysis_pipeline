function fitOutput = GU_stairsHistPlot(v1, varargin)
% ------------------input----------------------
% takes a vector input and generates a stairs histogram
% v1 is arranged in a cell array; MAX 5 arrays!!
% each cell should contain a vector
% -----------------------------------------------
% orange; blue; red; green; purple
% ------------------output----------------------
% fitOutput has three cell outputs for each vector
% 1. mu; 2. sigma; 3. A
% -----------------------------------------------
% overlaps can be visualized in matlab 2014a or below

% Gokul Upadhyayula, Feb 2015

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('v1');
ip.addParamValue('xlabel', '' , @isstr);
ip.addParamValue('ylabel', '' , @isstr);
ip.addParamValue('title', '' , @isstr);
ip.addParamValue('BinSize', 1 , @isnumeric);
ip.addParamValue('LineWidth', 1 , @isnumeric);
ip.addParamValue('FontSize', 22, @isnumeric);
ip.addParamValue('CreateFigure', true , @islogical);
ip.addParamValue('Fit', false , @islogical); % gaussian fit the data
ip.addParamValue('PercentileStart', 0, @isnumeric); % truncates the start data
ip.addParamValue('PercentileEnd', 99, @isnumeric); % truncates the end data
ip.addParamValue('ConstrainMeans', true , @islogical);
ip.addParamValue('ConstrainSD', true , @islogical);
ip.addParamValue('MinPercentile', 3, @isnumeric); % start search for mean for fitting
ip.addParamValue('ShowECDF', false, @islogical); % start search for mean for fitting
ip.parse(v1, varargin{:});

warning('off', 'all')
% warning
FS = ip.Results.FontSize;
maxPopulations = 10; % check 10 populations for the gaussian fitting
% hist colors
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO] = GU_colorCodes;
colorCell{1} = ceO;
colorCell{2} = ceB;
colorCell{3} = ceR;
colorCell{4} = ceG;
colorCell{5} = ceP;
colorCell{6} = cfO;
colorCell{7} = cfB;
colorCell{8} = cfR;
colorCell{9} = cfG;
colorCell{10} = cfP;
colorCell{11} = ceO;
colorCell{12} = cfO;

if ~iscell(v1)
    v1 = {v1};
end
% determine # of bins
vec1 = [];
for i = 1:numel(v1)
    if ~isempty(v1{i})
        if isrow(v1{i})
            vec1 = [vec1 v1{i}];
        else
            vec1 = [vec1 v1{i}'];
        end
    end
end
tmax = prctile(vec1,ip.Results.PercentileEnd)*1.1;
tmin = prctile(vec1,ip.Results.PercentileStart);%*0.9;
dt = ip.Results.BinSize;
t = tmin:dt:tmax;
nt = numel(t); % used for bins

if ip.Results.CreateFigure
    figure
end
hold on
for k = 1:numel(v1)
    if ~isempty(v1{k})
        if k <= 6
            yyaxis left
            set(gca, 'ycolor', 'k');
            histUT = hist(v1{k}, t);
            stairsXT(t, histUT/sum(histUT*dt), 'FaceColor', colorCell{k+5}, 'EdgeColor', colorCell{k}, 'LineWidth', ip.Results.LineWidth);
        else
            yyaxis left
            set(gca, 'ycolor', 'k');
            histUT = hist(v1{k}, t);
            stairsXT(t, histUT/sum(histUT*dt), 'FaceColor', colorCell{k}, 'EdgeColor', colorCell{k}, 'LineWidth', ip.Results.LineWidth);
        end
        % plot ecdf on the right axis
        if ip.Results.ShowECDF
            [f,x] = ecdf(v1{k});
            yyaxis right
            set(gca, 'ycolor', 'k');
            plot(x,f,':', 'Color',  colorCell{k}, 'LineWidth', ip.Results.LineWidth);
        end
        % fitting
        if ip.Results.Fit
            [SM_pdf, SM_pdf_x] =hist(v1{k}, nt); % generate centers and counts
            parfor n = 1:maxPopulations %checks 10 dist
                [~, ~, ~, ~, tBIC(n), ~] = fitGaussianMixture1D_v2(SM_pdf_x, SM_pdf,n,'ConstrainMeans', ip.Results.ConstrainMeans,'ConstrainSD',ip.Results.ConstrainSD, 'MinPercentile', ip.Results.MinPercentile,'Display', 'off');
            end
            [~,arg2] = min(tBIC);
            if arg2 == maxPopulations
                fprintf('Auto Minimization of BIC Failed... \n');
                fprintf('Please check the number of arguments tested and increase if appropriate\n');
            else
                for n = 1:arg2
                    [mu, sigma, A, ~, ~, ~] = fitGaussianMixture1D_v2(SM_pdf_x, SM_pdf,n,'ConstrainMeans', ip.Results.ConstrainMeans , 'ConstrainSD',ip.Results.ConstrainSD,'MinPercentile', ip.Results.MinPercentile,'Display', 'off');
                end
                
                % overlay fits
                SM_pdf_x = (min(v1{k})*0.8):(max(v1{k})/200):(max(v1{k})*1.20);
                sum_plot = zeros(1, numel(SM_pdf_x));
                
                for i = 1:1:arg2
                    normalized_plot=(normpdf(SM_pdf_x, mu(i), sigma(i)))*A(i);
                    sum_plot = sum_plot + normalized_plot;
                    hold on
                    plot(SM_pdf_x, normalized_plot, 'Color',colorCell{k+1}, 'LineWidth', ip.Results.LineWidth);
                end
                
                hold on
                plot (SM_pdf_x, sum_plot,'--','Color',colorCell{k+6}, 'LineWidth', ip.Results.LineWidth);
                
                fitOutput{k}{1} = mu;
                fitOutput{k}{2} = sigma;
                fitOutput{k}{3} = A;
            end
        else
            fitOutput{k} = [];
        end
        if ip.Results.CreateFigure
            yyaxis left
            set(gca,'fontsize',FS)
            xlabel(ip.Results.xlabel,'FontSize',FS, 'FontName', 'Helvetica');
            ylabel(ip.Results.ylabel,'FontSize',FS, 'FontName', 'Helvetica');
            title(ip.Results.title,'FontSize',FS, 'FontName', 'Helvetica');
            if ip.Results.ShowECDF
                hold on
                yyaxis right
                set(gca,'fontsize',FS)
                ylabel('Cumulative Frequency');
                ylim([0 1]);
                set(gca, 'Ytick', 0:0.2:1);
            end
        else
            yyaxis left
            xlabel(ip.Results.xlabel,'FontSize',FS, 'FontName', 'Helvetica');
            ylabel(ip.Results.ylabel,'FontSize',FS, 'FontName', 'Helvetica');
            title(ip.Results.title,'FontSize',FS, 'FontName', 'Helvetica');
            if ip.Results.ShowECDF
                hold on
                yyaxis right
                set(gca,'fontsize',FS)
                ylabel('Cumulative Frequency');
                ylim([0 1]);
                set(gca, 'Ytick', 0:0.2:1);
            end
        end
    end
end
