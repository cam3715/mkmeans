classdef mixedclust
    %MIXEDCLUST is a class for copmuting kmeans clustering
    %for data sets with numeric and categorical variables.
    %
    % Other m-files required:
    %     assignmentoptimal.m, Markus Buehren ...
    %         http://www.mathworks.com/matlabcentral/fileexchange/ ...
    %         6543-functions-for-the-rectangular-assignment-problem ...
    %         /content/assignmentoptimal.m
    %     This function implements the Hungarian Algorithm.
    %
    % Subfunctions:
    %     significance       (C) Ahmad Alsahaf - GNU GPL
    %     cluster_center     (C) Ahmad Alsahaf - GNU GPL
    %     algo_dist          (C) Ahmad Alsahaf - GNU GPL
    %     dist_to_center     (C) Ahmad Alsahaf - GNU GPL
    %
    %     These subfunctions can also be found as part of Ahmad Alsahaf's...
    %       amjams/mixedkmeans package on MATLAB Central/FileExchange.  They ...
    %       have been modified for use, and are inluded in mixedclust.m.
    %       http:///www.mathworks.com/matlabcentral/fileexchange
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
    
    
    
    %------------- BEGIN CODE --------------
    
    properties
        data
        m_idx
        k
        max_iter
        inputType
        trialsNo
        mixedClustering
        numericClustering
        origionalData
        significances
        data_discrete
        normalizedData
        tempvar
        all_dist
        silh_mean
        performance
        idx
    end
    
    methods
        function obj = mixedclust(data, k, max_iter, inputType,trialsNo)
            
            [dn, ~] = size(data);
            [~, ~] = size(inputType);
            if nargin < 4
                trialsNo = 1;
                inputType = [];
            elseif nargin < 3
                max_iter = 1000;
            elseif nargin < 2
                display('Not Enough Arguments')
            end
            obj.tempvar.dn = dn;
            obj.trialsNo = trialsNo;
            obj.data = data;
            obj.k = k;
            obj.max_iter = max_iter;
            obj.inputType = inputType;
            
            % replace NaN entrieies
            obj.data(isnan(obj.data)) = 1;
            obj = normalize(obj);
            obj = discretize(obj);
            obj = sigpairs(obj);
            obj = signif(obj);
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
        function obj = normalize(obj)
            obj.tempvar.m_distance = zeros(obj.tempvar.dn,obj.k,obj.trialsNo);
            obj.tempvar.n_distance = zeros(obj.tempvar.dn,obj.k,obj.trialsNo);
            
            
            
            %% Mixed KMeans
            obj.tempvar.silhouette_mixed_mean = zeros(1,obj.trialsNo);
            
            obj.m_idx = zeros(obj.tempvar.dn,obj.trialsNo);
            
            [obj.tempvar.n,obj.tempvar.m] = size(obj.data);
            obj.tempvar.idx_num = find(~obj.inputType);
            
            %% Normalize Numeric Data
            for i=1:numel(obj.tempvar.idx_num)
                obj.data(:,obj.tempvar.idx_num(i)) = (obj.data(:,obj.tempvar.idx_num(i)) - ...
                    repmat(min(obj.data(:,obj.tempvar.idx_num(i))),size(obj.data(:,obj.tempvar.idx_num(i))))) ...
                    /(max(obj.data(:,obj.tempvar.idx_num(i)))-min(obj.data(:,obj.tempvar.idx_num(i))));
            end
            obj.normalizedData = obj.data;
        end
        function obj = discretize(obj)
            %% Discretize Numeric Data
            obj.data_discrete = obj.data;
            %ensure k << N
            %             if obj.tempvar.dn >1000
            %                 max_k = round(obj.tempvar.dn/200);
            %             else
            max_k = 20;
            %             end
            obj.tempvar.idx_cat = find(-1*obj.inputType+1);
            for i=1:numel(obj.tempvar.idx_cat)
                silh_avg = zeros(max_k,1);
                data_num = obj.data(:,obj.tempvar.idx_cat(i));
                for k_iter=1:max_k
                    
                    [idx,~,~,D]=kmeans(data_num,k_iter+1,'dist', ...
                        'sqeuclidean','MaxIter',100,'Options',statset('UseParallel',1));
                    [Drow,~] = size(D);
                    silh = zeros(1,Drow);
                    for drow = 1:Drow
                        [a_drow,excl_D] = min(D(drow,:));
                        b_drow = min(D(drow,[1:(excl_D-1),(excl_D+1):end]));
                        silh(drow) = (b_drow-a_drow)/max(a_drow,b_drow);
                    end
                    %Ensure selection has >1 unique value
                    if numel(unique(idx))>1
                        silh_avg(k_iter) = mean(silh);
                    else
                        silh_avg(k_iter) = -10000000;
                    end
                    
                end
                [~,k_best] = max(silh_avg);
                k_best = k_best+1;
                obj.data_discrete(:,obj.tempvar.idx_cat(i)) = kmeans(obj.data(:,obj.tempvar.idx_cat(i)),k_best);
            end
        end
        function obj = sigpairs(obj)
            D = obj.data_discrete;
            
            % define the attribute, its unique values, and all unique pairs
            for i=1:obj.tempvar.m
                
                a = D(:,i);
                unique_a = find(accumarray(a+1,1))-1;
                all_pairs = nchoosek(unique_a,2);
                varname = ['var', num2str(i)];
                obj.tempvar.(varname).all_pairs = all_pairs;
                
            end
        end
        function obj = signif(obj)
            
            sigs = zeros(obj.tempvar.m,1);
            parfor i=1:obj.tempvar.m
                
                sigs(i) = significance(obj,i);
                
            end
            obj.significances = sigs;
        end
        function obj = mclust(obj)
            n = obj.tempvar.dn;
            for iMixed = 1:obj.trialsNo
                %try
                curr_idx = randi([1 obj.k],n,1);
                
                obj = algo_distance(obj);
                
                new_idx = zeros(n,1);
                
                count = 0;
                while(isequal(new_idx,curr_idx)==0 && count<obj.max_iter)
                    
                    if count>0
                        curr_idx = new_idx;
                    end
                    
                    all_centers = struct;
                    droppedACluster = 0;
                    for i=1:obj.k
                        curr_cluster = obj.data(curr_idx==i,:);
                        curr_center = cluster_center(curr_cluster,obj);
                        if curr_center.cluster_size < 1
                            droppedACluster = 1;
                        end
                        name = ['center_',sprintf('%03d',i)];
                        all_centers.(name) = curr_center;
                    end
                    
                    silh_c = zeros(1,n);
                    
                    if droppedACluster == 0
                        
                        for i=1:n
                            k_distances = zeros(obj.k,1);
                            data_i = obj.data(i,:);
                            for j=1:obj.k
                                name_now = ['center_',sprintf('%03d',j)];
                                center_now = all_centers.(name_now);
                                obj.tempvar.data_i = data_i;
                                obj.tempvar.center_now = center_now;
                                obj.tempvar.all_dist = obj.all_dist;
                                k_distances(j) = dist_to_center(obj);
                            end
                            
                            [~,new_idx(i)] = min(k_distances);
                            %                         m_distance(i,:,iMixed) = k_distances;
                            min1 = min(k_distances);
                            min2 = min(setdiff(k_distances(:),min(k_distances(:))));
                            silh_c(i) = (min2-min1)/max(min1,min2);
                        end
                    else
                        new_idx = randi([1 obj.k],n,1);
                    end
                    count = count+1;
                end
                
                idx = new_idx;
                obj.m_idx(:,iMixed) = idx;
                obj.silh_mean(iMixed) = mean(silh_c);
                %                 catch
                %                     fprintf('Error Non-existent field categorical. iMixed = %d', ...
                %                         iMixed);
                %                     display('- - - -  execution will continue  - - - -')
                %                     iMixed = iMixed-1;
                %                 end
            end
            
        end
        %% Significances
        function sig = significance(obj,idx)
            %SIGNIFICANCE: finds the significance of a categorical attribute or a
            %discretized version of a numerical attribute
            
            % input:
            %   D:   the dataset of all attributes
            %   idx: index of the attribute whose significance is to be found
            %
            % ouput:
            %   sig: the significance of the attribute
            %
            %
            % Copyright 2015 Ahmad Alsahaf
            % Research fellow, Politecnico di Milano
            % ahmadalsahaf@gmail.com
            
            % number of attributes
            D = obj.data_discrete;
            m = size(D,2);
            
            % define the attribute, its unique values, and all unique pairs
            a = D(:,idx);
            unique_a = find(accumarray(a+1,1))-1;
            varname = ['var', num2str(idx)];
            all_pairs = obj.tempvar.(varname).all_pairs;
            %Note nchoosek is impractical for n>15
            num_pairs = size(all_pairs,1);
            
            % the number of all delta distances
            num_delta = (m-1)*num_pairs;
            
            % find all deltas and average them
            feature_c = 1:m;  feature_c(idx)=[];   %complementary feature set
            delta_sum = 0;                         %initialize
            for i=1:num_pairs
                curr_pair = all_pairs(i,:);
                for j=1:(m-1)
                    % intialize distance
                    d = 0;
                    
                    % initalize support set
                    w = [];
                    w_c = [];
                    
                    % the number of categorical values in D(:,feature_c(j))
                    data_temp = D(:,feature_c(j))+1;
                    unique_j = find(accumarray(data_temp,1))-1;
                    vj = numel(unique_j);
                    
                    
                    % begin algorithm
                    for t = 1:vj
                        ut = unique_j(t);
                        
                        % locations
                        ut_in_aj = find(D(:,feature_c(j))==ut);
                        x_in_ai = find(a==curr_pair(1));
                        y_in_ai = find(a==curr_pair(2));
                        
                        % probabilities
                        p_ux = numel(x_in_ai(ismembc(x_in_ai,ut_in_aj))) ...
                            /numel(x_in_ai);
                        p_uy = numel(y_in_ai(ismembc(y_in_ai,ut_in_aj))) ...
                            /numel(y_in_ai);
                        
                        % conditions
                        if p_ux>= p_uy
                            w = [w;ut];          %update support set
                            d = d+p_ux;          %update distance
                        else
                            w_c = [w_c;p_uy];    %update complement support set
                            d = d+p_uy;          %update distance
                        end
                        
                        
                    end
                    delta = d-1;                      %restrict distance to [0,1]
                    delta_sum = delta_sum + delta;
                end
            end
            % find average distance, which is the significance
            sig = delta_sum/num_delta;
        end
        
        
        %% Algo Distance
        
        function obj = algo_distance(obj)
            % Copyright 2015 Ahmad Alsahaf
            % Research fellow, Politecnico di Milano
            % ahmadalsahaf@gmail.com
            data_discrete = obj.data_discrete;
            
            % data dimenionsality
            [~,m] = size(data_discrete);
            
            %intialize distance vector; which contains all distances between all pairs
            all_dist = [];
            
            for i = 1:m
                % define ai, the current attribute
                ai = data_discrete(:,i);
                
                % find all pairs of unique values in current feature
                unique_ai = find(accumarray(ai+1,1))-1;
                all_pairs = nchoosek(unique_ai,2);
                
                % find complement feature set
                feat_c = 1:m;  feat_c(i) = [];
                
                for j= 1:size(all_pairs,1)
                    % initialize sum and define current pair
                    sum_delta = 0;
                    curr_pair = all_pairs(j,:);
                    
                    % find distance between the pair for all Aj
                    for k = 1:m-1
                        % define aj
                        aj = data_discrete(:,feat_c(k));
                        
                        % update the sum
                        % intialize distance
                        d = 0;
                        
                        % initalize support set
                        w = [];
                        w_c = [];
                        
                        % the number of categorical values in aj
                        unique_j = find(accumarray(aj+1,1))-1;
                        
                        vj = numel(unique_j);
                        
                        % begin algorithm
                        for t = 1:vj
                            ut = unique_j(t);
                            
                            % locations
                            ut_in_aj = find(aj==ut);
                            x_in_ai = find(ai==curr_pair(1));
                            y_in_ai = find(ai==curr_pair(2));
                            
                            % probabilities
                            p_ux = numel(x_in_ai(ismembc(x_in_ai,ut_in_aj))) ...
                                /numel(x_in_ai);
                            p_uy = numel(y_in_ai(ismembc(y_in_ai,ut_in_aj))) ...
                                /numel(y_in_ai);
                            
                            % conditions
                            if p_ux>= p_uy
                                w = [w;ut];          %update support set
                                d = d+p_ux;          %update distance
                            else
                                w_c = [w_c;p_uy];    %update complement support set
                                d = d+p_uy;          %update distance
                            end
                            
                            
                        end
                        delta = d-1;                      %restrict distance to [0,1]
                        sum_delta = sum_delta + delta;
                    end
                    %         update the distance vector
                    sum_delta = sum_delta/(m-1);
                    
                    %         arranged as [attribute_idx,first_value(lower), ...
                    %         second_value(higher),distance];
                    pair_sorted = sort(curr_pair,'ascend');
                    all_dist = [all_dist; ...
                        i,pair_sorted(1),pair_sorted(2),sum_delta];
                    obj.all_dist = all_dist;
                end
            end
            
        end
        
        %% Cluster Center
        function [ center ] = cluster_center(cluster, obj)
            %CLUSTER_CENTER find cluster centers for mixed attributes
            
            % inputs:
            %     cluster:    the members of the cluster
            %     input_type: binary index indicating the type of attributes...
            %     (1 for categorical)
            %
            % output:
            %     center:     the center of the cluster
            % Copyright 2015 Ahmad Alsahaf
            % Research fellow, Politecnico di Milano
            % ahmadalsahaf@gmail.com
            
            % intialize a structure variable to save the centers
            
            input_type = obj.inputType;
            center = struct;
            
            % cluster dimensions, and numerical and categorical feature indices
            [n,~] = size(cluster);
            center.cluster_size = n;
            cat_idx = find(input_type);
            num_idx = find(~input_type);
            
            
            % find center for each numerical attribute
            for i=1:numel(num_idx)
                curr_att = cluster(:,num_idx(i));
                name = ['att_',sprintf('%03d',num_idx(i))];
                center.numerical.(name) = mean(curr_att);
            end
            
            % find center for each categorical attribute
            for i=1:numel(cat_idx)
                curr_att = cluster(:,cat_idx(i));
                name = ['att_',sprintf('%03d',cat_idx(i))];
                uniq_curr_att = find(accumarray(curr_att+1,1))-1;
                
                for j=1:numel(uniq_curr_att)
                    name_value = ['value_',sprintf('%03d',j)];
                    curr_value = uniq_curr_att(j);
                    count_value = numel(find(curr_att==curr_value));
                    center.categorical.(name).(name_value).value = curr_value;
                    center.categorical.(name).(name_value).count = count_value;
                end
            end
            
        end
        
        %% Distance to Center
        
        function theta = dist_to_center(obj)
            
            % dist_to_center: computes the the distance between a data point and a
            % cluster center
            
            % inputs:
            %     x:              a data point
            %     c:              a cluster center (structure)
            %     input_type:     binary index indicating attributes (1 = categorical)
            %     sig:            significance of all attributes in the dataset
            %     dist_all:       list of all distances of categorical values
            %
            % output:
            %     theta:  the distance between x and c
            %
            %
            % Copyright 2015 Ahmad Alsahaf
            % Research fellow, Politecnico di Milano
            % ahmadalsahaf@gmail.com
            
            x = obj.tempvar.data_i;
            c = obj.tempvar.center_now;
            input_type = obj.inputType;
            sig = obj.significances;
            dist_all = obj.tempvar.all_dist;
            
            % find indices
            cat_idx = find(input_type);
            num_idx = find(~input_type);
            
            % load cluster size
            cluster_size = c.cluster_size;
            
            % distance for numerical attributes
            
            % initialize numerical distance to zero
            sum_distance_numerical_v = zeros(1,numel(num_idx));
            
            % find distance for each numerical attribute and add to sum
            for i=1:numel(num_idx)
                d = x(num_idx(i));
                name = ['att_',sprintf('%03d',num_idx(i))];
                num_center = c.numerical.(name);
                curr_significance = sig(num_idx(i));
                curr_dist = (curr_significance*(d-num_center))^2;
                sum_distance_numerical_v(i) = curr_dist;
            end
            sum_distance_numerical = sum(sum_distance_numerical_v);
            
            % display(c)
            % initialize categorical distance to zero
            sum_distance_categorical = 0;
            
            % find distance for each categorical attribute and add to sum
            for i=1:numel(cat_idx)
                % access the current categorical attribute from structure
                name = ['att_',sprintf('%03d',cat_idx(i))];
                curr_att = c.categorical.(name);
                value_names = fieldnames(curr_att);
                
                % initialize sum for this categorical attribute
                sum_categorical_current = zeros(1,numel(value_names));
                
                % now access values within that attribute in the cluster
                
                for j=1:numel(value_names)
                    value_in_point = x(cat_idx(i));
                    value_in_cluster = curr_att.(value_names{j}).value;
                    count_in_cluster = curr_att.(value_names{j}).count;
                    
                    % find the distance from the list
                    sorted_values = sort([value_in_point,value_in_cluster],'ascend');
                    idx_dist = dist_all(:,1)==cat_idx(i)&dist_all(:,2) == ...
                        sorted_values(1) & dist_all(:,3) == sorted_values(2);
                    
                    % set distance to zero if value is equal to center, compute dist
                    % othewise (i.e. only update when different values
                    
                    if (sorted_values(1) ~= sorted_values(2))
                        sum_categorical_current(j) = ((1/cluster_size)*...
                            count_in_cluster*dist_all(idx_dist,4))^2;
                    end
                end
                sum_distance_categorical = sum(sum_categorical_current);
            end
            % overall distance
            theta = sum_distance_numerical + sum_distance_categorical;
        end
        
        
    end
end