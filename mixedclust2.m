classdef mixedclust
    %MIXEDCLUST is a class for comparing mixed and numeric kmeans clustering
    %for data sets with numeric and categorical variables.
    %
    %Inputs:
    %     data               - mxn matrix of n data points wiht m attributes
    %     k                  - number of means (cluster centers)
    %     max_iter           - maximum iterations for k-means
    %     inputType          - binary: attribute is (1= categorical 0= numeric)
    %     output             - given indicies for clusters in known structure
    %     trialsNo           - maximum number of means for clustering in the
    %                             discretation of numeric data
    %
    % Outputs:
    %   obj
    %         .data
    %         .k
    %         .max_iter
    %         .inputType
    %         .tiralsNo
    %         .mixedClustering
    %             .k
    %             .idx
    %             .silhouette
    %             .distance
    %             .trialsNo
    %         .numericClustering
    %             .k
    %             .idx
    %             .silhouette
    %             .distance
    %             .trialsNo
    %         .origionalData.data
    %         .data_discrete
    %         .normalizedData
    %         .tempvar
    %                 .dn
    %                 .m_distance
    %                 .n_distance
    %                 .waitlen
    %                 .silhouette_mixed_mean
    %                 .m_idx
    %                 .n
    %                 .idx_num
    %                 .data_i
    %                 .center_now
    %                 .all_dist
    %         .all_dist
    %
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
        k
        max_iter
        inputType
        trialsNo
        mixedClustering
        origionalData
        significances
        data_discrete
        normalizedData
        tempvar
        all_dist
    end
    
    methods
        function obj = mixedclust(data, k, max_iter, inputType,trialsNo)
            
            [dn, ~] = size(data);
            [~, ~] = size(inputType);
            if nargin < 4
                trialsNo = 2;
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
            obj = normalize(obj);
            obj = discretize(obj);
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
            
            
            obj.tempvar.waitlen = .05;
            
            %% Mixed KMeans
            obj.tempvar.silhouette_mixed_mean = zeros(1,obj.trialsNo);
            
            obj.tempvar.m_idx = zeros(obj.tempvar.dn,obj.trialsNo);
            
            tic

            
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
            max_k = 50;
            obj.tempvar.idx_cat = find(obj.inputType);
            for i=1:numel(obj.tempvar.idx_cat)
                silh_avg = zeros(max_k,1);
                data_num = obj.data(:,obj.tempvar.idx_cat(i));
                parfor k_iter=1:max_k
                    [~,~,~,D]=kmeans(data_num,k_iter+1,'dist', ...
                        'sqeuclidean','MaxIter',10000,'Options',statset('UseParallel',1));
                    [Drow,~] = size(D);
                    silh = zeros(1,Drow);
                    for drow = 1:Drow
                        [a_drow,excl_D] = min(D(drow,:));
                        b_drow = min(D(drow,[1:(excl_D-1),(excl_D+1):end]));
                        silh(drow) = (b_drow-a_drow)/max(a_drow,b_drow);
                    end
                    silh_avg(k_iter) = mean(silh);
                end
                [~,k_best] = max(silh_avg);
                k_best = k_best+1;
                obj.data_discrete(:,obj.tempvar.idx_cat(i)) = kmeans(obj.data(:,obj.tempvar.idx_cat(i)),k_best);
                if unique(obj.data_discrete(:,obj.tempvar.idx_cat(i))) < 2
                    display('not enough unique values')
                end
            end
        end
        function obj = cluster(obj)
            n = obj.tempvar.dn;
            iMixed = 1;
            while iMixed < obj.trialsNo+1
                %try
                curr_idx = randi([1 obj.k],n,1);
                
                data_discrete = obj.data_discrete;
                
                obj = algo_distance(obj);
                all_dist = obj.all_dist;
                
                new_idx = zeros(n,1);
                
                count = 0;
                while(isequal(new_idx,curr_idx)==0 && count<obj.max_iter)
                    
                    if count>0
                        curr_idx = new_idx;
                    end
                    
                    all_centers = struct;
                    for i=1:obj.k
                        curr_cluster = obj.data(curr_idx==i,:);
                        curr_center = cluster_center(curr_cluster,obj);
                        name = ['center_',sprintf('%03d',i)];
                        all_centers.(name) = curr_center;
                    end
                    
                    silh_c = zeros(1,n);
                    for i=1:n
                        k_distances = zeros(obj.k,1);
                        data_i = obj.data(i,:);
                        for j=1:obj.k
                            name_now = ['center_',sprintf('%03d',j)];
                            center_now = all_centers.(name_now);
                            obj.tempvar.data_i = data_i;
                            obj.tempvar.center_now = center_now;
                            obj.tempvar.all_dist = all_dist;
                            k_distances(j) = dist_to_center(obj);
                        end
                        
                        [~,new_idx(i)] = min(k_distances);
                        m_distance(i,:,iMixed) = k_distances;
                        min1 = min(k_distances);
                        min2 = min(setdiff(k_distances(:),min(k_distances(:))));
                        silh_c(i) = (min2-min1)/max(min1,min2);
                    end
                    count = count+1;
                end
                
                idx = new_idx;
                m_idx(:,iMixed) = idx;
                silhouette_mixed_mean(iMixed) = mean(silh_c);
                obj.tempvar.waitlen = obj.tempvar.waitlen + (.9-obj.tempvar.waitlen)/20;
                %                 catch
                %                     fprintf('Error Non-existent field categorical. iMixed = %d', ...
                %                         iMixed);
                %                     display('- - - -  execution will continue  - - - -')
                %                     iMixed = iMixed-1;
                %                 end
                iMixed = iMixed+1;
            end
            toc
        end
        %% Significances
        function sig = significance(D,idx)
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
            m = size(D,2);
            
            % define the attribute, its unique values, and all unique pairs
            a = D(:,idx);
            unique_a = find(accumarray(a+1,1))-1;
            all_pairs = nchoosek(unique_a,2);
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
                    unique_j = find(accumarray(D(:,feature_c(j))+1,1))-1;
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
                if unique(ai) < 2
                    display 'not enough unique values'
                    display(i);
                end
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
            sum_distance_numerical = 0;
            
            % find distance for each numerical attribute and add to sum
            for i=1:numel(num_idx)
                d = x(num_idx(i));
                name = ['att_',sprintf('%03d',num_idx(i))];
                num_center = c.numerical.(name);
                curr_significance = sig(num_idx(i));
                curr_dist = (curr_significance*(d-num_center))^2;
                sum_distance_numerical = sum_distance_numerical+curr_dist;
            end
            
            % display(c)
            % initialize categorical distance to zero
            sum_distance_categorical = 0;
            
            % find distance for each categorical attribute and add to sum
            
            for i=1:numel(cat_idx)
                % access the current categorical attribute from structure
                name = ['att_',sprintf('%03d',cat_idx(i))];
                curr_att = c.categorical.(name);
                
                
                % initialize sum for this categorical attribute
                sum_categorical_current = 0;
                
                % now access values within that attribute in the cluster
                value_names = fieldnames(curr_att);
                
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
                        sum_categorical_current = sum_categorical_current ...
                            + (1/cluster_size)*count_in_cluster*dist_all(idx_dist,4);
                    end
                end
                sum_distance_categorical = sum_distance_categorical ...
                    +(sum_categorical_current)^2;
            end
            
            % overall distance
            theta = sum_distance_numerical + sum_distance_categorical;
        end
        
        
    end
end