function [data_profile_means] = get_profile_means(data_struct,layer,n_length)

%% get fields in structure

fields = fieldnames(data_struct);

%% separate data into profiles

unique_dives = unique(data_struct.dive);
unique_profiles = 1:numel(unique_dives)*2;
for n = 1:numel(unique_dives)
    profile_index(n).vals = [unique_dives(n) unique_dives(n)];
end
profile_index = [profile_index.vals];

n_prof = 0;
for n = 1:numel(unique_profiles)
    
    if mod(n,2) == 1
        n_prof = n_prof+1;
        check = data_struct.dive == profile_index(n_prof) & data_struct.downup == 1 & data_struct.depth <= layer;
        for n_fields = 1:numel(fields)
            if eval(['numel(data_struct.',cell2mat(fields(n_fields)),')']) == n_length
                eval(['data_profile_means(n_prof).',cell2mat(fields(n_fields)),' = nanmean(data_struct.',cell2mat(fields(n_fields)),'(check));']);
            end
        end
    end
    if mod(n,2) == 0 
        n_prof = n_prof+1;        
        check = data_struct.dive == profile_index(n_prof) & data_struct.downup == 2 & data_struct.depth <= layer;
        for n_fields = 1:numel(fields)
            if eval(['numel(data_struct.',cell2mat(fields(n_fields)),')']) == n_length
                eval(['data_profile_means(n_prof).',cell2mat(fields(n_fields)),' = nanmean(data_struct.',cell2mat(fields(n_fields)),'(check));']);
            end
        end        
    end
end

%% same again but for top 10 m (considered surface)

n_prof = 0;
for n = 1:numel(unique_profiles)
    
    if mod(n,2) == 1
        n_prof = n_prof+1;
        check = data_struct.dive == profile_index(n_prof) & data_struct.downup == 1 & data_struct.depth <= 10;
        for n_fields = 1:numel(fields)
            if eval(['numel(data_struct.',cell2mat(fields(n_fields)),')']) == n_length
                eval(['data_profile_means(n_prof).',cell2mat(fields(n_fields)),'_surf = nanmean(data_struct.',cell2mat(fields(n_fields)),'(check));']);
            end
        end
    end
    if mod(n,2) == 0 
        n_prof = n_prof+1;        
        check = data_struct.dive == profile_index(n_prof) & data_struct.downup == 2 & data_struct.depth <= 10;
        for n_fields = 1:numel(fields)
            if eval(['numel(data_struct.',cell2mat(fields(n_fields)),')']) == n_length
                eval(['data_profile_means(n_prof).',cell2mat(fields(n_fields)),'_surf = nanmean(data_struct.',cell2mat(fields(n_fields)),'(check));']);
            end
        end        
    end
end

















end