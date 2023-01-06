function table_name = spmup_BIDS_QCtables(subjects, type)

% takes a subjects structure and create a table.tsv of QC metrics along
% with plotting the data
%
% FORMAT table_name = spmup_BIDS_QCtables(subjects, type)
%
% INPUT - subjects is a structure saved by run_spmup_bids
%         it is saved on disk as spmup_subjects_task-X.mat
%       - type is either 'anat' or 'fMRI'
%
% OUTPUT table_name is the name of the table.tsv
%             saved next to the subject structure
%
% Cyril Pernet
% ---------------------------
%  Copyright (C) SPMUP Team

supported = {'anat','fMRI'};
if all(~strcmpi(type,supported))
    error('unsupported type')
end

if ischar(subjects)
    filepath = fileparts(subjects);
    subjects = load(subjects);
    subjects = subjects.subjects;
else
    filepath = pwd;
end

switch lower(type)
    case 'anat'
        for s=size(subjects,2):-1:1
            Nsess(s) = length(subjects{s}.anat_qa);
            % sum(cellfun(@(x) isempty(x), subjects{s}.anat_qa));
        end

        for session=1:max(Nsess)
            % check it's not empty
            for s=size(subjects,2):-1:1
                N(s) = ~isempty(subjects{s}.anat_qa{session});
            end
            % extract metrics
            if any(N)
                AnatQA = table( cellfun(@(x) x.anat_qa{session}.SNR, subjects)',...
                                cellfun(@(x) x.anat_qa{session}.CNR, subjects)',...
                                cellfun(@(x) x.anat_qa{session}.FBER, subjects)',...
                                cellfun(@(x) x.anat_qa{session}.EFC, subjects)',...
                    'VariableNames',{'SNR','CNR','FBER','EFC'});
                AnatQA.Properties.Description = ['spmup AnatQC session ' num2str(session)];
                table_name{session} = [filepath filesep 'AnatQC_session-' num2str(session) '.tsv']; %#ok<*AGROW>
                writetable(AnatQA,table_name{session},...
                    'Delimiter','\t','FileType','text')
                spmup_plotqc(AnatQA,'new')
            end
        end

    case 'fmri'
        for s=size(subjects,2):-1:1
            Nsess(s) = length(subjects{s}.anat_qa);
            % sum(cellfun(@(x) isempty(x), subjects{s}.anat_qa));
        end

        for session=1:max(Nsess)
            table_index = 1;
            % check it's not empty
            for s=size(subjects,2):-1:1
                N(s) = ~isempty(subjects{s}.anat_qa{session});
            end
            % extract metrics
            if any(N)
                for s=size(subjects,2):-1:1
                    R(s) = size(subjects{s}.func_qa{session},1);
                end

                for run = 1:max(R)
                    fMRIQA = table( cellfun(@(x) x.func_qa{session}{run}.preproc_tSNR.GM, subjects)',...
                                    cellfun(@(x) x.func_qa{session}{run}.preproc_tSNR.WM, subjects)',...
                                    cellfun(@(x) x.func_qa{session}{run}.preproc_tSNR.CSF, subjects)',...
                                    cellfun(@(x) x.func_qa{session}{run}.preproc_tSNR.average, subjects)',...
                        'VariableNames',{   ['tSNR GM run ' num2str(run)],...
                                            ['tSNR WM run ' num2str(run)],...
                                            ['tSNR CSF run ' num2str(run)],...
                                            ['tSNR average run ' num2str(run)]});
                    if max(R) == 1
                        fMRIQA.Properties.Description = ['spmup tSNR session ' num2str(session)];
                        table_name{table_index} = [filepath filesep 'fMRIQC_session-' num2str(session) '.tsv'];
                    else
                        fMRIQA.Properties.Description = ['spmup tSNR session ' num2str(session) ' run ' num2str(run)];
                        table_name{table_index} = [filepath filesep 'fMRIQC_session-' num2str(session) '_run-' num2str(run) '.tsv'];
                    end
                    writetable(fMRIQA,table_name{table_index},'Delimiter','\t','FileType','text')
                    table_index = table_index+1;
                    spmup_plotqc(fMRIQA,'new');
                end
            end
        end
end

% clean-up empty cells
table_name(cellfun(@(x) isempty(x), table_name)) = [];

