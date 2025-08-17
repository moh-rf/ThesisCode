clear all; clc; close all;

fprintf('شروع پردازش...\n');

%%%%%%%%%%%%%%%%%%%%%%%% تنظیم مسیرها و پارامترها %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MyToolboxDir = fullfile('D:', 'payanname', 'Matlabcodes', 'bbci_toolbox'); 
WorkingDir   = fullfile('D:', 'payanname', 'Matlabcodes');              
NirsMyDataDir= fullfile('D:', 'NIRS');                                   
BehaviorDir  = fullfile('D:','Data_FNIRS','behavior','N-back','summary'); 

% بررسی وجود پوشه‌ها
if ~isfolder(MyToolboxDir), error('پوشه تولباکس BBCI یافت نشد: %s', MyToolboxDir); end
if ~isfolder(NirsMyDataDir), error('پوشه داده NIRS یافت نشد: %s', NirsMyDataDir); end
if ~isfolder(BehaviorDir), error('پوشه داده رفتاری یافت نشد: %s', BehaviorDir); end

cd(WorkingDir); 


fprintf('راه اندازی تولباکس BBCI...\n');
cd(MyToolboxDir);
try
    startup_bbci_toolbox('DataDir', NirsMyDataDir, 'TmpDir', tempdir, 'History', 0);
    addpath(genpath(pwd)); 
    cd(WorkingDir); 
    fprintf('تولباکس BBCI با موفقیت راه اندازی شد.\n');
catch ME
    cd(WorkingDir); 
    error('خطا در راه اندازی تولباکس BBCI: %s', ME.message);
end


% پارامترهای پردازش fNIRS
ival_epo  = [-5 60]*1000; % بازه زمانی اپوک حول مارکر شروع سری (میلی‌ثانیه)
ival_base = [-5 -2]*1000; % بازه زمانی بیس‌لاین (میلی‌ثانیه)
ival_hbo_extract = [10 20]; % بازه زمانی برای استخراج میانگین HbO (ثانیه)
markers_per_session = 9;   % فرض: 9 سری (مارکر) در هر جلسه
n_sessions = 3;

% لیست شرکت‌کنندگان
subdir_list_nirs = {'VP001-NIRS','VP002-NIRS','VP003-NIRS','VP004-NIRS','VP005-NIRS','VP006-NIRS','VP007-NIRS','VP008-NIRS','VP009-NIRS','VP010-NIRS','VP011-NIRS','VP012-NIRS','VP013-NIRS','VP014-NIRS','VP015-NIRS','VP016-NIRS','VP017-NIRS','VP018-NIRS','VP019-NIRS','VP020-NIRS','VP021-NIRS','VP022-NIRS','VP023-NIRS','VP024-NIRS','VP025-NIRS','VP026-NIRS'};
subdir_list_behav = {'VP001-EEG','VP002-EEG','VP003-EEG','VP004-EEG','VP005-EEG','VP006-EEG','VP007-EEG','VP008-EEG','VP009-EEG','VP010-EEG','VP011-EEG','VP012-EEG','VP013-EEG','VP014-EEG','VP015-EEG','VP016-EEG','VP017-EEG','VP018-EEG','VP019-EEG','VP020-EEG','VP021-EEG','VP022-EEG','VP023-EEG','VP024-EEG','VP025-EEG','VP026-EEG'};
num_participants = length(subdir_list_nirs);

% لیست‌های موقت برای ذخیره نتایج fNIRS همه شرکت‌کنندگان
AllSubs_SubjectList_fNIRS = {};
AllSubs_SessionList_fNIRS = [];
AllSubs_ChannelList_fNIRS = {};
AllSubs_ConditionList_fNIRS = {};
AllSubs_HbOList_fNIRS = [];

% ماتریس‌های ذخیره نتایج رفتاری (Subject, Condition: 1=0back, 2=2back, 3=3back, Session)
all_accuracy_sess = nan(num_participants, 3, n_sessions);
all_react_sess = nan(num_participants, 3, n_sessions);
% *** ایجاد متغیر جدید برای ذخیره انحراف معیار RT ***
all_react_sd_sess = nan(num_participants, 3, n_sessions);

%%%%%%%%%%%%%%%%%%%%%%%% حلقه اصلی پردازش شرکت‌کنندگان %%%%%%%%%%%%%%%%%%%%%%
for vp = 1:num_participants
    fprintf('\n===== پردازش شرکت‌کننده %d: %s =====\n', vp, subdir_list_nirs{vp});

    % --- بخش پردازش fNIRS ---
    fprintf('--- پردازش fNIRS ---\n');
    loadDir_nirs = fullfile(NirsMyDataDir, subdir_list_nirs{vp});
    if ~isfolder(loadDir_nirs), warning('پوشه NIRS برای %s یافت نشد، پرش.', subdir_list_nirs{vp}); continue; end
    cd(loadDir_nirs);

    try
        load cnt_nback.mat;
        load mrk_nback.mat;
        load mnt_nback.mat;
        fprintf('فایل‌های cnt, mrk, mnt برای n-back بارگذاری شد.\n');
    catch ME
        warning('خطا در بارگذاری فایل‌های n-back برای %s: %s. پرش از fNIRS.', subdir_list_nirs{vp}, ME.message);
        cd(WorkingDir);
        continue; % پرش به شرکت‌کننده بعدی اگر فایل‌ها ناقص باشند
    end

    % صحت‌سنجی ساختار مارکر
    expected_markers = markers_per_session * n_sessions;
    if ~isfield(mrk_nback,'y') || ~ismatrix(mrk_nback.y) || size(mrk_nback.y, 2) ~= expected_markers || size(mrk_nback.y, 1) ~= 3
        warning('VP %d: ساختار mrk_nback.y (%dx%d) با انتظار 3x%d همخوانی ندارد! پرش از fNIRS.', vp, size(mrk_nback.y,1), size(mrk_nback.y,2), expected_markers);
        cd(WorkingDir);
        continue;
    end
    if ~isfield(mrk_nback,'className') || ~iscell(mrk_nback.className) || length(mrk_nback.className) ~= 3
        warning('VP %d: ساختار mrk_nback.className با انتظار 1x3 cell همخوانی ندارد! پرش از fNIRS.', vp);
        cd(WorkingDir);
        continue;
    end

    % فیلتر کردن داده پیوسته
    fprintf('اعمال فیلتر Butterworth...\n');
    [b,a] = butter(6, 0.2 / (cnt_nback.oxy.fs/2), 'low');
    cnt_nback.oxy = proc_filtfilt(cnt_nback.oxy, b, a);
    cnt_nback.deoxy = proc_filtfilt(cnt_nback.deoxy, b, a);

    % لیست‌های موقت برای نتایج fNIRS *این* شرکت‌کننده
    SubjectList_fNIRS_vp = {};
    SessionList_fNIRS_vp = [];
    ChannelList_fNIRS_vp = {};
    ConditionList_fNIRS_vp = {};
    HbOList_fNIRS_vp = [];

    for sess = 1:n_sessions % حلقه برای هر Session
        fprintf('  پردازش fNIRS - جلسه %d...\n', sess);
        % انتخاب مارکرهای این Session
        start_idx = (sess - 1) * markers_per_session + 1;
        end_idx = sess * markers_per_session;
        indices_sess = start_idx:end_idx;

        mrk_this_session = struct();
        mrk_this_session.time = mrk_nback.time(indices_sess);
        mrk_this_session.y = mrk_nback.y(:, indices_sess);
        mrk_this_session.className = mrk_nback.className;
        mrk_this_session.fs = cnt_nback.oxy.fs; % استفاده از fs داده cnt

        % اپوک بندی فقط برای این session
        epo_sess.oxy = proc_segmentation(cnt_nback.oxy, mrk_this_session, ival_epo);
        epo_sess.deoxy = proc_segmentation(cnt_nback.deoxy, mrk_this_session, ival_epo);

        % تصحیح بیس‌لاین
        epo_sess.oxy = proc_baseline(epo_sess.oxy, ival_base);
        epo_sess.deoxy = proc_baseline(epo_sess.deoxy, ival_base);

        % میانگین گیری بر اساس کلاس در این session
        epo_avg_sess.oxy = proc_average(epo_sess.oxy, 'Stats', 1);
        epo_avg_sess.deoxy = proc_average(epo_sess.deoxy, 'Stats', 1);

        % استخراج HbO از اپوک‌های میانگین
        t = epo_avg_sess.oxy.t / 1000;
        idx_range = find(t >= ival_hbo_extract(1) & t <= ival_hbo_extract(2));
        channels = epo_avg_sess.oxy.clab;
        condition_names_in_avg = epo_avg_sess.oxy.className;
        num_classes_in_avg = size(epo_avg_sess.oxy.x, 3);

        if isempty(idx_range)
            warning('VP %d, Session %d: بازه زمانی %d-%d ثانیه در اپوک یافت نشد.', vp, sess, ival_hbo_extract(1), ival_hbo_extract(2));
            continue; % پرش از این session اگر بازه زمانی معتبر نباشد
        end

        if length(condition_names_in_avg) ~= num_classes_in_avg
            warning('VP %d, Session %d: عدم تطابق نام کلاس (%d) و بعد داده (%d)!', vp, sess, length(condition_names_in_avg), num_classes_in_avg);
            continue;
        end

        for cond_idx = 1:num_classes_in_avg
            current_condition_name = condition_names_in_avg{cond_idx};
            condition_label = 'Unknown'; % پیش‌فرض
            if contains(current_condition_name, '0')
                condition_label = '0-back';
            elseif contains(current_condition_name, '2')
                condition_label = '2-back';
            elseif contains(current_condition_name, '3')
                condition_label = '3-back';
            else
                warning('VP %d, Session %d: نام کلاس "%s" تشخیص داده نشد.', vp, sess, current_condition_name);
            end

            for ch = 1:length(channels)
                mean_val = mean(epo_avg_sess.oxy.x(idx_range, ch, cond_idx));

                SubjectList_fNIRS_vp{end+1,1} = subdir_list_nirs{vp}(1:6); % e.g., "VP001"
                ChannelList_fNIRS_vp{end+1,1} = channels{ch};
                ConditionList_fNIRS_vp{end+1,1} = condition_label;
                SessionList_fNIRS_vp(end+1,1) = sess;
                HbOList_fNIRS_vp(end+1,1) = mean_val;
            end
        end
        fprintf('  جلسه %d fNIRS تمام شد.\n', sess);
    end % پایان حلقه session fNIRS

    % اضافه کردن نتایج این شرکت‌کننده به لیست کلی
    AllSubs_SubjectList_fNIRS = [AllSubs_SubjectList_fNIRS; SubjectList_fNIRS_vp];
    AllSubs_SessionList_fNIRS = [AllSubs_SessionList_fNIRS; SessionList_fNIRS_vp];
    AllSubs_ChannelList_fNIRS = [AllSubs_ChannelList_fNIRS; ChannelList_fNIRS_vp];
    AllSubs_ConditionList_fNIRS = [AllSubs_ConditionList_fNIRS; ConditionList_fNIRS_vp];
    AllSubs_HbOList_fNIRS = [AllSubs_HbOList_fNIRS; HbOList_fNIRS_vp];
    fprintf('--- پردازش fNIRS برای %s تمام شد ---\n', subdir_list_nirs{vp});

    % --- بخش پردازش رفتاری ---
    fprintf('--- پردازش رفتاری ---\n');
    vpDir_behav = fullfile(BehaviorDir, subdir_list_behav{vp});
    if ~isfolder(vpDir_behav), warning('پوشه رفتاری برای %s یافت نشد، پرش.', subdir_list_behav{vp}); cd(WorkingDir); continue; end
    cd(vpDir_behav);

    try
        load summary1.mat; load summary2.mat; load summary3.mat;
        sessions_data = {summary1, summary2, summary3};
        fprintf('فایل‌های summary1, 2, 3 بارگذاری شد.\n');
    catch ME
        warning('خطا در بارگذاری فایل‌های summary برای %s: %s. پرش از پردازش رفتاری.', subdir_list_behav{vp}, ME.message);
        cd(WorkingDir);
        continue;
    end

    for s = 1:n_sessions % حلقه برای هر Session رفتاری
        fprintf('  پردازش رفتاری - جلسه %d...\n', s);
        curr_summary = sessions_data{s};

        % --- کد محاسبه flag (کپی شده از کد شما با کمی اصلاح) ---
        if ~isfield(curr_summary, 'flag') || ~isfield(curr_summary, 'nback')
            warning('VP %d, Session %d: فیلدهای flag یا nback در summary وجود ندارند.', vp, s);
            continue; % پرش از این session
        end
        curr_summary.flag_org = curr_summary.flag; % نگه داشتن نسخه اصلی
        curr_summary.flag = zeros(size(curr_summary.flag_org)); % مقداردهی اولیه flag هدف
        target_0back = 1; target_2back = 2; target_3back = 3;

        for i = 1:size(curr_summary.flag,1) % loop over blocks?
            current_nback_type = curr_summary.nback(i);
            switch current_nback_type
                case 0
                    curr_summary.flag(i,:) = target_0back; % All trials are potential targets
                case 2
                    for j = 3:size(curr_summary.flag,2) % loop over trials within block
                        if curr_summary.flag_org(i, j-2) == curr_summary.flag_org(i, j)
                            curr_summary.flag(i,j) = target_2back; % Target trial
                            % else: flag remains 0 (non-target)
                        end
                    end
                case 3
                    for j = 4:size(curr_summary.flag,2)
                        if curr_summary.flag_org(i, j-3) == curr_summary.flag_org(i, j)
                            curr_summary.flag(i,j) = target_3back; % Target trial
                            % else: flag remains 0 (non-target)
                        end
                    end
                otherwise
                    warning('VP %d, Session %d: نوع n-back ناشناخته (%d).', vp, s, current_nback_type);
            end
        end
        % --- پایان محاسبه flag ---

        % --- محاسبه Accuracy و RT و SDRT برای این session ---
        results_all_trials = curr_summary.result(:); % تبدیل به بردار ستونی
        flags_all_trials = curr_summary.flag(:);     % تبدیل به بردار ستونی
        rt_all_trials = curr_summary.reaction_time(:); % تبدیل به بردار ستونی
        resp_all_trials = curr_summary.response(:);   % تبدیل به بردار ستونی

        conditions_to_calc = {'0-back', '2-back', '3-back'};
        targets = [target_0back, target_2back, target_3back];
        condition_indices = [1, 2, 3]; % 1 for 0-back, 2 for 2-back, 3 for 3-back

        for c_idx = 1:length(conditions_to_calc)
            target_val = targets(c_idx);
            cond_map_idx = condition_indices(c_idx);

            % انتخاب trial های مربوط به این condition (بر اساس flag هدف)
            relevant_trials_idx = (flags_all_trials == target_val);

            % محاسبه Accuracy (بدون تغییر)
            if any(relevant_trials_idx)
                correct_responses = sum(results_all_trials(relevant_trials_idx) == 1);
                total_relevant_trials = sum(relevant_trials_idx);
                accuracy = correct_responses / total_relevant_trials;
            else
                accuracy = NaN; % هیچ trial مرتبطی یافت نشد
            end
            all_accuracy_sess(vp, cond_map_idx, s) = accuracy;

            % محاسبه Mean RT و STD RT برای پاسخ‌های صحیح
            if strcmp(conditions_to_calc{c_idx}, '0-back')
                correct_rt_idx = relevant_trials_idx & (resp_all_trials==7 | resp_all_trials==8) & (results_all_trials == 1);
            else % 2-back and 3-back
                correct_rt_idx = relevant_trials_idx & (resp_all_trials==7) & (results_all_trials == 1);
            end

            % استخراج زمان واکنش‌های معتبر
            valid_rts = rt_all_trials(correct_rt_idx);
            % حذف مقادیر NaN احتمالی از RTها
            valid_rts = valid_rts(~isnan(valid_rts));

            if ~isempty(valid_rts)
                mean_rt = mean(valid_rts);
                % *** محاسبه انحراف معیار ***
                if length(valid_rts) > 1
                    sd_rt = std(valid_rts);
                else
                    sd_rt = 0; % یا NaN
                end
            else
                mean_rt = NaN;
                sd_rt = NaN;
            end
            all_react_sess(vp, cond_map_idx, s) = mean_rt;
            % *** ذخیره انحراف معیار ***
            all_react_sd_sess(vp, cond_map_idx, s) = sd_rt;

            % نمایش نتایج session
            fprintf('      %s: Accuracy = %.2f%%, Mean RT = %.3f ms, SD RT = %.3f ms\n', ...
                conditions_to_calc{c_idx}, accuracy * 100, mean_rt, sd_rt);
        end
        fprintf('  جلسه %d رفتاری تمام شد.\n', s);
    end % پایان حلقه session رفتاری
    fprintf('--- پردازش رفتاری برای %s تمام شد ---\n', subdir_list_behav{vp});

    cd(WorkingDir); % بازگشت به پوشه کاری اصلی
end % ===== پایان حلقه شرکت‌کنندگان =====

fprintf('\n===== پردازش تمام شرکت‌کنندگان تمام شد. شروع ادغام نهایی... =====\n');

%%%%%%%%%%%%%%%%%%%%%%%% ادغام نهایی و ساخت جدول %%%%%%%%%%%%%%%%%%%%%%%%%%

% ساخت جدول خام fNIRS از لیست‌های کلی
HbOTable_raw = table(AllSubs_SubjectList_fNIRS, AllSubs_SessionList_fNIRS, AllSubs_ChannelList_fNIRS, AllSubs_ConditionList_fNIRS, AllSubs_HbOList_fNIRS, ...
    'VariableNames', {'Subject_str','Session','Channel','Condition','Mean_HbO'});

if isempty(HbOTable_raw)
    error('جدول fNIRS خام خالی است. خطایی در پردازش رخ داده است.');
end

% تبدیل Subject_str به شماره فرد (1 تا 26)
subjectNums_fNIRS = cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), HbOTable_raw.Subject_str);

% ایجاد جدول نهایی
HbOTable_final = HbOTable_raw;
HbOTable_final.Subject = subjectNums_fNIRS; % اضافه کردن ستون عددی Subject
HbOTable_final = removevars(HbOTable_final, 'Subject_str'); % حذف ستون رشته‌ای

% داده‌های دموگرافیک
AgeList = [26, 26, 25, 30, 26, 25, 32, 24, 27, 22, ...
    27, 33, 27, 22, 23, 30, 17, 25, 25, 26, ...
    28, 22, 28, 31, 29, 22];
GenderList = {'female', 'female', 'female', 'male', 'female', ...
    'female', 'female', 'female', 'female', 'male', ...
    'male', 'male', 'male', 'female', 'male', 'male', ...
    'female', 'female', 'female', 'male', 'male', ...
    'female', 'female', 'female', 'female', 'female'};

% محاسبه تعداد ردیف‌ها برای هر فرد در جدول نهایی
rows_per_subject_map = containers.Map('KeyType','double','ValueType','double');
unique_subjects_in_table = unique(HbOTable_final.Subject); % گرفتن افراد موجود در جدول نهایی
for i = 1:length(unique_subjects_in_table)
     subj_num = unique_subjects_in_table(i);
     rows_per_subject_map(subj_num) = sum(HbOTable_final.Subject == subj_num);
end


% ساخت لیست‌های توسعه‌یافته Age و Gender
ExpandedAgeList = [];
ExpandedGenderList = {};
sorted_subj_nums = sort(unique_subjects_in_table); % استفاده از افراد موجود

for i = 1:length(sorted_subj_nums)
    subj_num = sorted_subj_nums(i);
    if subj_num >= 1 && subj_num <= length(AgeList) && isKey(rows_per_subject_map, subj_num)
        num_rows = rows_per_subject_map(subj_num);
        ExpandedAgeList = [ExpandedAgeList; repelem(AgeList(subj_num), num_rows, 1)];
        ExpandedGenderList = [ExpandedGenderList; repelem(GenderList(subj_num), num_rows, 1)];
    else
        warning('شماره شرکت‌کننده %d نامعتبر است یا در جدول fNIRS یافت نشد.', subj_num);
         if isKey(rows_per_subject_map, subj_num)
             num_rows = rows_per_subject_map(subj_num);
             ExpandedAgeList = [ExpandedAgeList; repelem(NaN, num_rows, 1)];
             ExpandedGenderList = [ExpandedGenderList; repelem({''}, num_rows, 1)];
         end
    end
end

% اضافه کردن ستون‌های Age و Gender
if length(ExpandedAgeList) == height(HbOTable_final)
    HbOTable_final.Age = ExpandedAgeList;
    HbOTable_final.Gender = ExpandedGenderList;
else
     warning('طول لیست‌های Age/Gender با جدول نهایی مطابقت ندارد. بررسی کنید.');
     % در صورت عدم تطابق، ممکن است نیاز به روش دیگری برای اضافه کردن باشد
     HbOTable_final.Age = nan(height(HbOTable_final),1); % اضافه کردن ستون خالی
     HbOTable_final.Gender = repmat({''}, height(HbOTable_final),1);
end


% اضافه کردن ستون‌های Accuracy و MeanRT و SDRT ***
nRowsFinal = height(HbOTable_final);
AccuracyList = nan(nRowsFinal, 1);
RTList = nan(nRowsFinal, 1);
% *** اضافه شده: لیست برای SDRT ***
SDRTList = nan(nRowsFinal, 1);

condition_map = containers.Map({'0-back', '2-back', '3-back'}, {1, 2, 3});

for i = 1:nRowsFinal
    subj_num = HbOTable_final.Subject(i);
    sess_num = HbOTable_final.Session(i);
    cond_str = HbOTable_final.Condition{i}; % نام قبلی Condition

    % *** تغییر نام ستون Condition به n_back در اینجا انجام می‌شود ***
    if isKey(condition_map, cond_str)
        cond_idx = condition_map(cond_str);

        if subj_num >= 1 && subj_num <= size(all_accuracy_sess, 1) && ...
                cond_idx >= 1 && cond_idx <= size(all_accuracy_sess, 2) && ...
                sess_num >= 1 && sess_num <= size(all_accuracy_sess, 3)

            AccuracyList(i) = all_accuracy_sess(subj_num, cond_idx, sess_num) * 100; % درصد
            RTList(i) = all_react_sess(subj_num, cond_idx, sess_num); % فرض: میلی‌ثانیه
            % *** اضافه شده: خواندن SDRT ***
            SDRTList(i) = all_react_sd_sess(subj_num, cond_idx, sess_num); % فرض: میلی‌ثانیه
        else
            warning('اندیس نامعتبر هنگام دسترسی به داده رفتاری: Subject=%d, CondIdx=%d, Session=%d (ردیف %d)', subj_num, cond_idx, sess_num, i);
        end
    else
        warning('Condition "%s" در ردیف %d برای مپینگ رفتاری یافت نشد.', cond_str, i);
    end
end

HbOTable_final.Accuracy = AccuracyList;
HbOTable_final.MeanRT = RTList;
% *** اضافه شده: اضافه کردن ستون SDRT ***
HbOTable_final.SDRT = SDRTList;

% تغییر نام ستون‌ها مطابق مقاله هدف (Cakar & Yavuz) و نیازهای نمودار R
fprintf('تغییر نام ستون‌ها...\n');
HbOTable_final.Properties.VariableNames{'Channel'}   = 'Indices';
HbOTable_final.Properties.VariableNames{'Condition'} = 'n_back'; % <<< نام Condition به n_back تغییر کرد
HbOTable_final.Properties.VariableNames{'Mean_HbO'}  = 'value';

% مرتب‌سازی ستون‌ها (اختیاری - SDRT اضافه شده)
desired_order = {'Subject', 'Session', 'Age', 'Gender', 'n_back', 'Indices', 'value', 'Accuracy', 'MeanRT', 'SDRT'}; % SDRT اضافه شد
if all(ismember(desired_order, HbOTable_final.Properties.VariableNames))
    HbOTable_final = HbOTable_final(:, desired_order);
else
    warning('برخی ستون‌های مورد نظر برای مرتب‌سازی در جدول نهایی وجود ندارند.');
    disp('ستون‌های موجود:');
    disp(HbOTable_final.Properties.VariableNames);
end


%%%%%%%%%%%%%%%%%%%%%%%% نمایش و ذخیره نهایی %%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n===== جدول نهایی ادغام شده (شامل SDRT) =====\n');
disp(head(HbOTable_final, 10)); % نمایش 10 ردیف اول

output_filename = 'HbO_Behavior.csv'; % نام فایل جدید
output_fullpath = fullfile(WorkingDir, output_filename);
try
    writetable(HbOTable_final, output_fullpath);
    fprintf('\nجدول نهایی با موفقیت در فایل زیر ذخیره شد:\n%s\n', output_fullpath);
catch ME
    fprintf('\nخطا در ذخیره جدول نهایی:\n%s\n', ME.message);
end

fprintf('\n===== پردازش کامل شد =====\n');

% --- پایان ---