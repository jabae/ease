%% setup path
addpath(genpath(fullfile(EASE_dir, 'packages', 'microns_phase1_nda')));
addpath(genpath(fullfile(EASE_dir, 'packages', 'pipeline', 'matlab')));
addpath(genpath(fullfile(EASE_dir, 'packages', 'ta3')));
addpath(genpath(fullfile(EASE_dir, 'packages', 'ta3p100')));

%% database user account
database_list = {'ninai.cluster-chjk7zcxhsgn.us-east-1.rds.amazonaws.com:3306', ...
    '127.0.0.1:3306'};

% select database
if ~exist('dj_host', 'var')
    fprintf('************ SELECT A DATABASE ************\n');
    for m=1:length(database_list)
        fprintf('%d: %s\n', m, database_list{m});
    end
    fprintf('\n');
    temp = input('database ID: ');
    dj_host = database_list{temp};
    fprintf('\n**************** Done ********************\n');
end
fprintf('You are going to connect to a database\n\t%s.\n', dj_host); 

% type user names
fprintf('Now type your login information\n', dj_host);

if ~exist('dj_username', 'var')
    dj_username = input('username: ', 's');
end
if ~exist('dj_password', 'var')
    dj_password = input('password: ', 's');
end

%connect to the database
try
    setenv('DJ_HOST', dj_host)
    setenv('DJ_USER', dj_username)
    setenv('DJ_PASS', dj_password)
    dj.conn()
    fprintf('Database connected\n'); 
    clear dj_host dj_username dj_password; 
catch
    error('No connection with the selected database.');
end
