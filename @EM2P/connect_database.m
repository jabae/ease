function [] = connect_database(obj, databases, dj_username, dj_password)
%% connect to a database for fetching data
%{
%}

%% inputs
%{
	obj: type; description
	databases: cell array ; each element is one database host (string) 
	dj_username: string; username
	dj_password: string; user password 
%}

%% outputs
%{
%}

%% Author
%{
	Pengcheng Zhou 
	Columbia Unviersity, 2019
	zhoupc2018@gmail.com
	GPL-3.0 License 
%}



if obj.dj_connected
    fprintf('the database has been connected. \nif you want to reconnect, set ease.dj_connected=false first.\n');
    return;
end

%% database user account
if ~exist('databases_list', 'var')
    databases = obj.databases_list; 
else
    databases = union(obj.databases_list, databases); 
end

% select database
if ~exist('dj_host', 'var')
    fprintf('************ SELECT A DATABASE ************\n');
    for m=1:length(databases)
        fprintf('%d: %s\n', m, databases{m});
    end
    fprintf('\n');
    temp = input('database ID: ');
    dj_host = databases{temp};
    fprintf('\n**************** Done ********************\n');
end
fprintf('You are going to connect to a database\n\t%s.\n', dj_host);

% type user names
fprintf('Now type your login information\n');

if ~exist('dj_username', 'var')
    dj_username = input('username: ', 's');
end
if ~exist('dj_password', 'var')
    dj_password = input('password: ', 's');
    clc; 
end

%connect to the database
try
    setenv('DJ_HOST', dj_host)
    setenv('DJ_USER', dj_username)
    setenv('DJ_PASS', dj_password)
    dj.conn()
    fprintf('Database connected\n');
    clear dj_host dj_username dj_password;
    obj.dj_connected = true;
catch
    error('No connection with the selected database.');
end
