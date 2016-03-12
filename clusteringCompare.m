classdef clusteringCompare
    %CLUSTERINGCOMPARE Summary of this class goes here
    %   Detailed explanation goes here
    %
    % Author: Camden Glenn Bock
    % 598 Bates College, Lewistion, ME 04240
    % cbock@bates.edu, camdenbock@gmail.com
    % http://www.camdenbock.com
    % December 2015; Last Revision: 12/30/2015
    %
    %% Copyright (C) 2016  Camden Bock - GPL v. 3.0
    %
    % 'This program is liscensed under GPL v3.0'
    % 'This program is modified from Ahmad Alsahaf`s package: ...
    % amjams/mixedkmeans'.
    %
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <http://www.gnu.org/licenses/>.
    %
    %
    % 'This program comes with ABSOLUTELY NO WARRANTY;' ...
    % 'for details type view source. This is free software, and' ...
    % 'you are welcome to redistribute it display under certain' ...
    % 'conditions; see <http://www.gnu.org/licenses/>', ...
    % ('Copyright (C) 2016  Camden Bock, Bates College'))
    
    
    
    %% ------------- BEGIN CODE --------------
    
    properties
        mixedClust
        numericClust
        output
        trialsNo
        inputType
        data
        tempvar
    end
    methods
        function obj = clusteringCompare(filename, catAttributes, ncols, ...
                startRow, trialsNo)
            if nargin < 5
                trialsNo = 1;
            elseif nargin < 4
                startRow = 1;
            elseif nargin < 3
                display('not enough inputs')
            end
            obj.trialsNo = trialsNo;
            obj = dataimport(obj,filename, catAttributes, ncols, startRow);
            [n,m] = size(obj.data);
            obj.mixedClust = mixedclust(obj.data, obj.tempvar.k, 1000, obj.inputType,trialsNo);
            obj.mixedClust = mclust(obj.mixedClust);
            obj.mixedClust.idx = obj.mixedClust.m_idx(:,:,1);
            obj.numericClust = struct;
            obj.numericClust.idx = zeros(m,obj.trialsNo);
            obj.numericClust.D = zeros(m,obj.tempvar.k,obj.trialsNo);
            obj.numericClust.silh = zeros(n,obj.trialsNo);
            obj.numericClust.avg_silh = zeros(1,obj.trialsNo);
            for i=1:obj.trialsNo
                [idx,~,~,D] = kmeans(obj.data,obj.tempvar.k);
                obj.numericClust.idx = idx;
                obj.numericClust.dist(:,:,i) = D;
                silh = zeros(n,1);
                for j = 1:n
                    x = sort(D(j,:));
                    silh(j) = (x(2)-x(1))/x(2);
                end
                obj.numericClust.silh(:,i) = silh;
                obj.numericClust.avg_silh(i) = mean(silh);
            end
            obj = compareOut(obj);
        end
        function obj = compareOut(obj)
            obj.numericClust.performance = zeros(1,obj.trialsNo);
            obj.mixedClust.performance = zeros(1,obj.trialsNo);
            for nameidx=1:2
                if nameidx==1
                    name = 'mixedClust';
                elseif nameidx==2
                    name = 'numericClusst';
                end
                for i=1:obj.trialsNo
                    idx = obj.(name).idx(:,:,i);
                    k = numel(unique(obj.output));
                    ErrorMatrix = zeros(k);
                    output_values = unique(obj.output);
                    for emCol = 1:k
                        for emRow = 1:k
                            output_emRow = output_values(emRow);
                            for oRow = 1:length(obj.output)
                                if (obj.output(oRow) ~= output_emRow) && (idx(oRow) == emCol);
                                    ErrorMatrix(emCol,emRow) ...
                                        = ErrorMatrix(emCol,emRow)+1;
                                end
                            end
                        end
                    end
                    [mEM, nEM] = size(ErrorMatrix);
                    if mEM~=nEM
                        display('Warning, matrix must be square');
                    end
                    [~, count] = assignmentoptimal(ErrorMatrix);
                    obj.(name).performance(i) = 1-count/length(idx);
                end
            end
        end
        function obj = dataimport(obj,filename, catAttributes, ncols, startRow)
            
            delimiter = ',';
            endRow = inf;
            
            %% Format string for each line of text:  Autmoated for user input
            inputBlock = ('%f');
            formatSpec = char(1:(2*ncols + 8));
            for i = 2:2:(2*ncols)
                formatSpec((i-1):i) = inputBlock;
            end
            formatSpec((2*ncols+1):(2*ncols+8))=('%[^\n\r]');
            
            %% Open the text file.
            fileID = fopen(filename,'r');
            
            %% Read columns of data according to format string.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, ...
                'Delimiter', delimiter, 'HeaderLines', startRow(1)-1, ...
                'ReturnOnError', false);
            for block=2:length(startRow)
                frewind(fileID);
                dataArrayBlock = textscan(fileID, formatSpec, ...
                    endRow(block)-startRow(block)+1, 'Delimiter', delimiter, ...
                    'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
                for col=1:length(dataArray)
                    dataArray{col} = [dataArray{col};dataArrayBlock{col}];
                end
            end
            
            %% Close the text file.
            fclose(fileID);
            
            %% Create output variable
            datasource = [dataArray{1:end-1}];
            
            %% Options for alsahaf_mixed_kmeans
            obj.data = datasource(:,1:(end-1));
            obj.output = datasource(:,end);
            if min(min(obj.output))==0
                obj.output = obj.output + 1;
            elseif min(min(obj.output))~=1
                disp('Error: Datasource may not have proper categorical assignments')
            end
            
            obj.tempvar.k = length(unique(obj.output));
            
            [~,dc] = size(obj.data);
            obj.inputType = zeros(1,dc);
            for q=1:length(catAttributes)
                obj.inputType(catAttributes(q)) = 1;
                if min(min(obj.data(:,catAttributes(q))))== 0
                    obj.data(:,catAttributes(q)) = obj.data(:,catAttributes(q)) + 1;
                elseif min(min(obj.data(:,catAttributes(q))))~=1
                    disp('Error: Datasource may not have proper categorical assignments')
                end
            end
            obj.tempvar.inputType = obj.inputType;
        end
        function visulaizeNum(obj, trialnum)
            if nargin < 2
                trialnum = 1;
            end
            pointClusterVis(obj.numericClust, trialnum)
        end
        function visualizeMix(obj, trialnum)
            if nargin < 2
                trialnum = 1;
            end
            pointClusterVis(obj.mixedClust, trialnum)
        end
    end
    
end