DROP TABLE IF EXISTS users_auth;
DROP TABLE IF EXISTS users_meta;
DROP TABLE IF EXISTS users_affiliation;
DROP TABLE IF EXISTS projects;
DROP TABLE IF EXISTS projects_roles;
DROP TABLE IF EXISTS data;
DROP TABLE IF EXISTS data_deletion;
DROP TABLE IF EXISTS data_table_gene;
DROP TABLE IF EXISTS results;
DROP TABLE IF EXISTS analysis;
DROP TABLE IF EXISTS analysis_author;
DROP TABLE IF EXISTS analysis_inputs;
DROP TABLE IF EXISTS universities;
DROP TABLE IF EXISTS university_domains;

CREATE TABLE users_auth (
    id_user TEXT UNIQUE NOT NULL PRIMARY KEY,  -- id; backend
    name_user TEXT UNIQUE NOT NULL,  -- username; visible to user
    password TEXT NOT NULL,  -- this contains hash method, salt, and pass hash
    email TEXT UNIQUE NOT NULL,
    name TEXT,  -- human name, "John Smith"
    orcid TEXT UNIQUE,  -- linked ORCID
    creation_time_user TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    ip TEXT,
    confirmed BIT DEFAULT 0,  -- whether user has confirmed their email address
    confirmed_datetime TIMESTAMP
);

CREATE TABLE users_affiliation (
    id_user TEXT NOT NULL,
    affiliation TEXT NOT NULL,
    CONSTRAINT affiliation_key PRIMARY KEY (id_user, affiliation),
    FOREIGN KEY (id_user) REFERENCES users(id_user)
);

CREATE TABLE projects (
    id_project TEXT UNIQUE NOT NULL,  -- id for this project
    id_user TEXT NOT NULL,  -- owner of the project
    name_project TEXT NOT NULL,  -- title of project, input by user
    description TEXT,  -- description of project, input by user
    num_files INT DEFAULT 0,
    num_members INT DEFAULT 1,
    creation_time_project TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (id_user) REFERENCES users_auth(id_user)
);

CREATE TABLE projects_roles (
    id_project TEXT UNIQUE NOT NULL,
    id_user TEXT NOT NULL,
    role TEXT NOT NULL DEFAULT "member",
    data_read BIT DEFAULT 0,
    data_write BIT DEFAULT 0,
    FOREIGN KEY (id_project) REFERENCES projects(id_project),
    CONSTRAINT proj_key PRIMARY KEY (id_project, id_user)
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
    is_published BIT DEFAULT 0,  -- whether this data has been published; prevents auto deletion
    FOREIGN KEY (id_user) REFERENCES users(id_user),
    FOREIGN KEY (id_project) REFERENCES projects(id_project),
    FOREIGN KEY (id_results) REFERENCES results(id_results)
);

CREATE TABLE data_deletion (
    id_data TEXT UNIQUE NOT NULL PRIMARY KEY,
    id_user TEXT NOT NULL,  -- owner/uploader
    id_project TEXT NOT NULL,  -- project that this data is associated with
    file_path TEXT NOT NULL,  -- path of the data on server
    data_type TEXT NOT NULL,
    file_hash TEXT NOT NULL,  -- sha3-256 hash of the data
    id_results TEXT DEFAULT NULL,  -- pointer to results, if any
    data_uploaded_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    id_user_deleter TEXT NOT NULL,  -- user who marked data for deletion,
    data_deletion_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    deletion_path TEXT NOT NULL,  -- path where data is staged for deletion
    FOREIGN KEY (id_data) REFERENCES data(id_data),
    FOREIGN KEY (id_project) REFERENCES projects(id_project)
);

CREATE TABLE results(
    id_result TEXT UNIQUE NOT NULL PRIMARY KEY,  -- unique ID for each file produced by analysis
    id_project TEXT NOT NULL,
    id_analysis TEXT,  -- ID specific to analysis
    file_path TEXT NOT NULL,  -- where the result can be found
    result_type TEXT NOT NULL,
    description TEXT,  -- description of result
    title TEXT,
    is_interactive BIT NOT NULL DEFAULT 0,
    FOREIGN KEY (id_project) REFERENCES projects(id_project),
    FOREIGN KEY (id_analysis) REFERENCES analysis(id_analysis)
);

CREATE TABLE analysis(
    id_analysis TEXT UNIQUE NOT NULL PRIMARY KEY,  -- unique id for each pipeline
    id_project TEXT NOT NULL,  -- project with which this pipeline is associated
    analysis_name TEXT UNIQUE NOT NULL,  -- name used for display
    node_file TEXT NOT NULL,  -- file containing node configuration
    parameter_file TEXT NOT NULL,  -- file containing pipeline object
    analysis_hash NOT NULL,  -- hash of the node_file + parameter_file, for checking if analysis has changed since validation
    is_validated BIT NOT NULL DEFAULT 0,
    initial_pscs_version TEXT NOT NULL DEFAULT 0,
    FOREIGN KEY (id_project) REFERENCES projects(id_project)
);

CREATE TABLE analysis_author(
    id_analysis TEXT UNIQUE NOT NULL,
    id_user TEXT NOT NULL,
    FOREIGN KEY (id_analysis) REFERENCES analysis(id_analysis),
    FOREIGN KEY (id_user) REFERENCES users(id_user),
    CONSTRAINT author_key PRIMARY KEY (id_user, id_analysis)
);

CREATE TABLE analysis_inputs(
    id_input TEXT UNIQUE NOT NULL PRIMARY KEY,  -- unique id for each input - not used
    id_analysis TEXT NOT NULL,  -- analysis with which the input is associated
    node_id TEXT NOT NULL,  -- id of the node within the pipeline, not within DB
    node_name TEXT NOT NULL,  -- displayed name of the node
    FOREIGN KEY (id_analysis) REFERENCES analysis(id_analysis)
);

CREATE TABLE universities(
  id_university TEXT UNIQUE NOT NULL PRIMARY KEY,
  name TEXT,
  alpha_two_code TEXT,  -- country code
  state_province TEXT,
  country TEXT,
  web_page TEXT
);

CREATE TABLE university_domains(
  id_university TEXT NOT NULL,
  university_domain TEXT NOT NULL,
  FOREIGN KEY (id_university) REFERENCES universities(id_university),
  CONSTRAINT university_key PRIMARY KEY (id_university, university_domain)
);