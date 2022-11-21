DROP TABLE IF EXISTS users;
DROP TABLE IF EXISTS projects;
DROP TABLE IF EXISTS data;
DROP TABLE IF EXISTS data_table_gene;
DROP TABLE IF EXISTS results;
DROP TABLE IF EXISTS analysis;

CREATE TABLE users (
    id_user TEXT UNIQUE NOT NULL PRIMARY KEY,
    name_user TEXT UNIQUE NOT NULL,
    password TEXT NOT NULL,
    creation_time_user TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE projects (
    id_project TEXT UNIQUE NOT NULL,  -- id for this project
    id_user TEXT NOT NULL,
    name_project TEXT NOT NULL,
    description TEXT,
    role TEXT NOT NULL DEFAULT "member",
    num_files INT DEFAULT 0,
    num_members INT DEFAULT 1,
    creation_time_project TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT proj_key PRIMARY KEY (id_project, id_user),
    FOREIGN KEY (id_user) REFERENCES users(id_user)
);

CREATE TABLE data (
    id_data TEXT UNIQUE NOT NULL PRIMARY KEY,
    id_user TEXT NOT NULL,  -- owner/uploader
    id_project TEXT NOT NULL,  -- project that this data is associated with
    file_path TEXT NOT NULL,  -- path of the data
    data_type TEXT NOT NULL,
    file_hash TEXT NOT NULL,  -- sha3-256 hash of the data
    id_results TEXT DEFAULT NULL,  -- pointer to results, if any
    data_uploaded_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (id_user) REFERENCES users(id_user),
    FOREIGN KEY (id_project) REFERENCES projects(id_project),
    FOREIGN KEY (id_results) REFERENCES results(id_results)
);

CREATE TABLE results(
    id_result TEXT UNIQUE NOT NULL PRIMARY KEY,  -- unique ID for each file produced by analysis
    id_project TEXT NOT NULL,
    id_analysis TEXT,  -- ID specific to analysis
    file_path TEXT NOT NULL,  -- where the result can be found
    result_type TEXT NOT NULL,
    description TEXT,  -- description of result
    title TEXT,
    is_interactive BIT NOT NULL DEFAULT 0
);

CREATE TABLE analysis(
    id_analysis TEXT UNIQUE NOT NULL PRIMARY KEY,  -- unique id for each pipeline
    id_user TEXT NOT NULL,  -- creator of pipeline
    id_project TEXT NOT NULL,  -- project with which this pipeline is associated
    parameter_file TEXT NOT NULL  -- file containing pipeline object
);