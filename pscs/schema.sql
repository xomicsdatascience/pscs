PRAGMA foreign_keys = ON;

CREATE TABLE users_auth (
    id_user TEXT UNIQUE NOT NULL PRIMARY KEY,  -- id; backend
    name_user TEXT UNIQUE NOT NULL,  -- username; visible to user
    password TEXT NOT NULL,  -- this contains hash method, salt, and pass hash
    email TEXT UNIQUE NOT NULL,
    confirmed_phi BOOL DEFAULT 0,  -- unless otherwise specified, user has not accepted phi terms
    confirmed_datause BOOL DEFAULT 0,  -- unless otherwise specified, user has not accepted datause terms
    name TEXT,  -- human name, "John Smith"
    orcid TEXT UNIQUE,  -- linked ORCID
    creation_time_user TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    ip TEXT,
    confirmed BIT DEFAULT 0,  -- whether user has confirmed their email address
    confirmed_datetime TIMESTAMP
);

CREATE TABLE users_confirmation (  -- for keeping track of confirmation emails that have been sent out
    id_user TEXT UNIQUE NOT NULL PRIMARY KEY,  -- id of user
    num_sent INT DEFAULT 1,  -- number of times the confirmation email has been sent
    last_sent TIMESTAMP DEFAULT CURRENT_TIMESTAMP,  -- last time the email was sent
    FOREIGN KEY (id_user) REFERENCES users_auth(id_user) ON DELETE CASCADE
);  -- don't bother keeping records of this table

CREATE TABLE users_affiliation (
    id_user TEXT NOT NULL,
    affiliation TEXT NOT NULL,
    affiliation_order INTEGER,
    CONSTRAINT affiliation_key PRIMARY KEY (id_user, affiliation),
    FOREIGN KEY (id_user) REFERENCES users_auth(id_user) ON DELETE CASCADE
);

CREATE TABLE admin_user (
    id_user TEXT NOT NULL,  -- admin users;
    FOREIGN KEY (id_user) REFERENCES users_auth(id_user) ON DELETE CASCADE
);

CREATE TABLE projects (
    id_project TEXT UNIQUE NOT NULL,  -- id for this project
    id_user TEXT NOT NULL,  -- owner of the project
    name_project TEXT NOT NULL,  -- title of project, input by user
    description TEXT,  -- description of project, input by user
    num_files INT DEFAULT 0,
    num_members INT DEFAULT 1,
    creation_time_project TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    is_published BIT DEFAULT 0,  -- whether project is public
    is_peer_review BIT DEFAULT 0,  -- whether project is under peer review
    is_on_hold BIT DEFAULT 0,  -- whether project is transition to public/peer, awaiting user confirmation
    date_published TIMESTAMP,
    FOREIGN KEY (id_user) REFERENCES users_auth(id_user) ON DELETE CASCADE
);

CREATE TABLE projects_deletion (
    id_project TEXT NOT NULL,  -- id for this project
    id_user TEXT NOT NULL,  -- owner of the project
    name_project TEXT NOT NULL,  -- title of project, input by user
    description TEXT,  -- description of project, input by user
    num_files INT DEFAULT 0,
    num_members INT DEFAULT 1,
    creation_time_project TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    is_published BIT DEFAULT 0,
    is_peer_review BIT DEFAULT 0,
    is_on_hold BIT DEFAULT 0,
    date_published TIMESTAMP,
    deletion_time_projects TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

-- NOTE: This would be better done as a trigger that selects all columns automatically instead of listing them here.
CREATE TRIGGER stage_project_deletion
  BEFORE DELETE ON projects
  FOR EACH ROW
  BEGIN
    INSERT INTO projects_deletion(id_project, id_user, name_project, description, num_files, num_members, creation_time_project, is_published, is_peer_review, is_on_hold, date_published)
    VALUES(OLD.id_project, OLD.id_user, OLD.name_project, OLD.description, OLD.num_files, OLD.num_members, OLD.creation_time_project, OLD.is_published, OLD.is_peer_review, OLD.is_on_hold, OLD.date_published);
  END;

CREATE TABLE projects_roles (
    id_project TEXT,
    id_user TEXT NOT NULL,
    role TEXT NOT NULL DEFAULT 'member',
    data_read BIT DEFAULT 0,  -- whether user can read data uploaded to project
    data_write BIT DEFAULT 0,  -- whether user can upload new data
    analysis_read BIT DEFAULT 0,  -- whether user can load pipeline in editor
    analysis_write BIT DEFAULT 0,  -- whether user can add new pipelines to project
    analysis_execute BIT DEFAULT 0,  -- whether user can run pipelines on project data
    project_management BIT DEFAULT 0,  -- whether user can add users
    CONSTRAINT proj_key PRIMARY KEY (id_project, id_user),
    FOREIGN KEY (id_project) REFERENCES projects(id_project) ON DELETE CASCADE,
    FOREIGN KEY (id_user) REFERENCES users_auth(id_user) ON DELETE CASCADE
);

CREATE TABLE projects_roles_deletion (
    id_project TEXT,
    id_user TEXT NOT NULL,
    role TEXT NOT NULL DEFAULT 'member',
    data_read BIT DEFAULT 0,  -- whether user can read data uploaded to project
    data_write BIT DEFAULT 0,  -- whether user can upload new data
    analysis_read BIT DEFAULT 0,  -- whether user can load pipeline in editor
    analysis_write BIT DEFAULT 0,  -- whether user can add new pipelines to project
    analysis_execute BIT DEFAULT 0,  -- whether user can run pipelines on project data
    project_management BIT DEFAULT 0,  -- whether user can add users
    deletion_time_projects_roles TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT proj_key PRIMARY KEY (id_project, id_user)
);

CREATE TRIGGER stage_project_roles_deletion
  BEFORE DELETE ON projects_roles
  FOR EACH ROW
  BEGIN
    INSERT INTO projects_roles_deletion(id_project, id_user, role, data_read, data_write, analysis_read, analysis_write, analysis_execute, project_management)
    VALUES(OLD.id_project, OLD.id_user, OLD.role, OLD.data_read, OLD.data_write, OLD.analysis_read, OLD.analysis_write, OLD.analysis_execute, OLD.project_management);
  END;

CREATE TABLE data (
    id_data TEXT UNIQUE NOT NULL PRIMARY KEY,
    id_user TEXT NOT NULL,  -- owner/uploader
    id_project TEXT NOT NULL,  -- project that this data is associated with
    file_path TEXT NOT NULL,  -- path of the data
    data_type TEXT NOT NULL,
    file_hash TEXT NOT NULL,  -- sha3-256 hash of the data
    data_uploaded_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    is_published BIT DEFAULT 0,  -- whether this data has been published; prevents auto deletion
    is_peer_review BIT DEFAULT 0,  -- whether the data is under peer review
    file_name VARCHAR(100) NOT NULL DEFAULT 'data',
    FOREIGN KEY (id_user) REFERENCES users_auth(id_user) ON DELETE CASCADE,
    FOREIGN KEY (id_project) REFERENCES projects(id_project) ON DELETE CASCADE
);

CREATE TABLE data_deletion (
  id_data TEXT NOT NULL PRIMARY KEY,
  id_user TEXT NOT NULL,  -- owner/uploader
  id_project TEXT NOT NULL,  -- project that this data is associated with
  file_path TEXT NOT NULL,  -- path of the data
  data_type TEXT NOT NULL,
  file_hash TEXT NOT NULL,  -- sha3-256 hash of the data
  data_uploaded_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  is_published BIT DEFAULT 0,  -- whether this data has been published; prevents auto deletion
  is_peer_review BIT DEFAULT 0,  -- whether the data is under peer review
  deletion_time_data TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TRIGGER stage_data_deletion
  BEFORE DELETE ON data
  FOR EACH ROW
  BEGIN
    INSERT INTO data_deletion(id_data, id_user, id_project, file_path, data_type, data_type, file_hash, data_uploaded_time, is_published, is_peer_review)
    VALUES(OLD.id_data, OLD.id_user, OLD.id_project, OLD.file_path, OLD.data_type, OLD.data_type, OLD.file_hash, OLD.data_uploaded_time, OLD.is_published, OLD.is_peer_review);
  END;

CREATE TABLE results(
    id_result TEXT UNIQUE NOT NULL PRIMARY KEY,  -- unique ID for each file produced by analysis
    id_project TEXT NOT NULL,
    id_analysis TEXT,  -- ID specific to analysis
    file_path TEXT NOT NULL,  -- where the result can be found
    file_name VARCHAR(200),  -- display-friendly name
    result_type TEXT NOT NULL,
    description TEXT,  -- description of result
    title TEXT,
    is_interactive BIT NOT NULL DEFAULT 0,
    interactive_tag VARCHAR(3) DEFAULT NULL,
    is_published BIT DEFAULT 0,
    is_peer_review BIT DEFAULT 0,
    FOREIGN KEY (id_project) REFERENCES projects(id_project) ON DELETE CASCADE,
    FOREIGN KEY (id_analysis) REFERENCES analysis(id_analysis) ON DELETE CASCADE
);

CREATE TABLE results_deletion(
    id_result TEXT NOT NULL PRIMARY KEY,  -- unique ID for each file produced by analysis
    id_project TEXT NOT NULL,
    id_analysis TEXT,  -- ID specific to analysis
    file_path TEXT NOT NULL,  -- where the result can be found
    result_type TEXT NOT NULL,
    description TEXT,  -- description of result
    title TEXT,
    is_interactive BIT NOT NULL DEFAULT 0,
    is_published BIT DEFAULT 0,
    is_peer_review BIT DEFAULT 0,
    deletion_time_results TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);


CREATE TRIGGER stage_results_deletion
  BEFORE DELETE ON results
  FOR EACH ROW
  BEGIN
    INSERT INTO results_deletion(id_result, id_project, id_analysis, file_path, result_type, description, title, is_interactive, is_published, is_peer_review)
    VALUES(OLD.id_result, OLD.id_project, OLD.id_analysis, OLD.file_path, OLD.result_type, OLD.description, OLD.title, OLD.is_interactive, OLD.is_published, OLD.is_peer_review);
  END;

CREATE TABLE analysis(
    id_analysis TEXT UNIQUE NOT NULL PRIMARY KEY,  -- unique id for each pipeline
    id_project TEXT NOT NULL,  -- project with which this pipeline is associated
    analysis_name TEXT NOT NULL,  -- name used for display
    node_file TEXT NOT NULL,  -- file containing node configuration
    parameter_file TEXT NOT NULL,  -- file containing pipeline object
    analysis_hash NOT NULL,  -- hash of the node_file + parameter_file, for checking if analysis has changed since validation
    is_validated BIT NOT NULL DEFAULT 0,
    initial_pscs_version TEXT NOT NULL DEFAULT 0,
    is_published BIT DEFAULT 0,
    is_peer_review BIT DEFAULT 0,
    FOREIGN KEY (id_project) REFERENCES projects(id_project) ON DELETE CASCADE
);

CREATE TABLE analysis_deletion(
  id_analysis TEXT NOT NULL PRIMARY KEY,  -- unique id for each pipeline
  id_project TEXT NOT NULL,  -- project with which this pipeline is associated
  analysis_name TEXT NOT NULL,  -- name used for display
  node_file TEXT NOT NULL,  -- file containing node configuration
  parameter_file TEXT NOT NULL,  -- file containing pipeline object
  analysis_hash NOT NULL,  -- hash of the node_file + parameter_file, for checking if analysis has changed since validation
  is_validated BIT NOT NULL DEFAULT 0,
  initial_pscs_version TEXT NOT NULL DEFAULT 0,
  is_published BIT DEFAULT 0,
  is_peer_review BIT DEFAULT 0,
  deletion_time_analysis TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TRIGGER stage_analysis_deletion
  BEFORE DELETE ON analysis
  FOR EACH ROW
  BEGIN
    INSERT INTO analysis_deletion(id_analysis, id_project, analysis_name, node_file, parameter_file, analysis_hash, is_validated, initial_pscs_version, is_peer_review, is_published)
    VALUES(OLD.id_analysis, OLD.id_project, OLD.analysis_name, OLD.node_file, OLD.parameter_file, OLD.analysis_hash, OLD.is_validated, OLD.initial_pscs_version, OLD.is_peer_review, OLD.is_published);
  END;

CREATE TABLE analysis_author(
    id_analysis TEXT UNIQUE NOT NULL,
    id_user TEXT NOT NULL,
    FOREIGN KEY (id_analysis) REFERENCES analysis(id_analysis) ON DELETE CASCADE,
    FOREIGN KEY (id_user) REFERENCES users_auth(id_user) ON DELETE CASCADE,
    CONSTRAINT author_key PRIMARY KEY (id_user, id_analysis)
);

CREATE TABLE analysis_author_deletion(
    id_analysis TEXT NOT NULL,
    id_user TEXT NOT NULL,
    deletion_time_analysis_author TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TRIGGER stage_analysis_author_deletion
  BEFORE DELETE ON analysis_author
  FOR EACH ROW
  BEGIN
    INSERT INTO analysis_author_deletion(id_analysis, id_user)
    VALUES(OLD.id_analysis, OLD.id_user);
  END;

CREATE TABLE analysis_inputs(  -- specifies and describes which nodes of an analysis are input
    id_input TEXT UNIQUE NOT NULL PRIMARY KEY,  -- unique id for each input - not used
    id_analysis TEXT NOT NULL,  -- analysis with which the input is associated
    node_id TEXT NOT NULL,  -- id of the node within the pipeline, not within DB
    node_name TEXT NOT NULL,  -- displayed name of the node
    FOREIGN KEY (id_analysis) REFERENCES analysis(id_analysis) ON DELETE CASCADE
);

CREATE TABLE analysis_interactive_tags(
    interactive_tag VARCHAR(3) NOT NULL
);

CREATE TABLE analysis_inputs_deletion(  -- specifies and describes which nodes of an analysis are input
    id_input TEXT UNIQUE NOT NULL PRIMARY KEY,  -- unique id for each input - not used
    id_analysis TEXT NOT NULL,  -- analysis with which the input is associated
    node_id TEXT NOT NULL,  -- id of the node within the pipeline, not within DB
    node_name TEXT NOT NULL,  -- displayed name of the node
    deletion_time_analysis_inputs TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TRIGGER stage_analysis_inputs_deletion
  BEFORE DELETE ON analysis_author
  FOR EACH ROW
  BEGIN
    INSERT INTO analysis_author_deletion(id_analysis, id_user)
    VALUES(OLD.id_analysis, OLD.id_user);
  END;

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
  FOREIGN KEY (id_university) REFERENCES universities(id_university) ON DELETE CASCADE,
  CONSTRAINT university_key PRIMARY KEY (id_university, university_domain)
);

CREATE TABLE posts(
  id_post TEXT UNIQUE NOT NULL PRIMARY KEY,
  id_user TEXT NOT NULL,
  title TEXT NOT NULL,
  body TEXT NOT NULL,
  date_created TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  FOREIGN KEY (id_user) REFERENCES users_auth(id_user)
);

CREATE TABLE submitted_jobs(  -- logging submitted jobs
  id_job TEXT UNIQUE NOT NULL PRIMARY KEY,  -- unique id for job
  submitted_resource TEXT NOT NULL,  -- computing resource to which the job was submitted
  resource_job_id TEXT NOT NULL,  -- job id on the resource
  id_user TEXT NOT NULL,  -- submitter
  id_project TEXT NOT NULL,  -- project for which this was submitted
  id_analysis TEXT NOT NULL,  -- analysis that was submitted
  server_response TEXT,  -- text response from the server
  remote_results_directory TEXT NOT NULL, -- directory where results are stored
  date_submitted TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
  stdout_log VARCHAR(1024) DEFAULT NULL, -- stdout log of the result
  stderr_log VARCHAR(1024) DEFAULT NULL, -- stderr log of the result
  is_complete BIT DEFAULT 0,
  date_completed TIMESTAMP
);

CREATE TABLE submitted_jobs_deletion(  -- logging submitted jobs
   id_job TEXT NOT NULL PRIMARY KEY,  -- unique id for job
   submitted_resource TEXT NOT NULL,  -- computing resource to which the job was submitted
   resource_job_id TEXT NOT NULL,  -- job id on the resource
   id_user TEXT NOT NULL,  -- submitter
   id_project TEXT NOT NULL,  -- project for which this was submitted
   id_analysis TEXT NOT NULL,  -- analysis that was submitted
   server_response TEXT,  -- text response from the server
   remote_results_directory TEXT NOT NULL, -- directory where results are stored
   date_submitted TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
   is_complete BIT DEFAULT 0,
   date_completed TIMESTAMP,
   deletion_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TRIGGER stage_submitted_jobs_deletion
  BEFORE DELETE ON submitted_jobs
  FOR EACH ROW
  BEGIN
    INSERT INTO submitted_jobs_deletion(id_job, submitted_resource, resource_job_id, id_user, id_project, id_analysis, server_response, remote_results_directory, date_submitted, is_complete, date_completed)
    VALUES(OLD.id_job, OLD.submitted_resource, OLD.resource_job_id, OLD.id_user, OLD.id_project, OLD.id_analysis, OLD.server_response, OLD.remote_results_directory, OLD.date_submitted, OLD.is_complete, OLD.date_completed);
  END;

CREATE TABLE submitted_data(  -- for knowing which data was submitted with a job
  id_job TEXT NOT NULL,
  id_data TEXT NOT NULL,
  node_name TEXT NOT NULL,  -- node to which the data is submitted
  FOREIGN KEY (id_job) REFERENCES submitted_jobs(id_job) ON DELETE CASCADE,
  FOREIGN KEY (id_data) REFERENCES data(id_data) ON DELETE CASCADE,
  CONSTRAINT submitted_data_key PRIMARY KEY (id_job, id_data, node_name)
);

CREATE TABLE submitted_data_deletion(  -- for knowing which data was submitted with a job
   id_job TEXT NOT NULL,
   id_data TEXT NOT NULL,
   node_name TEXT NOT NULL,  -- node to which the data is submitted
   deletion_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
   CONSTRAINT submitted_data_key PRIMARY KEY (id_job, id_data, node_name)
);

CREATE TRIGGER stage_submitted_data_deletion
    BEFORE DELETE ON submitted_data
    FOR EACH ROW
    BEGIN
        INSERT INTO submitted_data_deletion(id_job, id_data, node_name)
        VALUES(OLD.id_job, OLD.id_data, OLD.node_name);
    END;

CREATE TABLE projects_peer_review(
    id_project TEXT NOT NULL,
    peer_password TEXT NOT NULL,  -- password supplied to journal/reviewers
    FOREIGN KEY (id_project) REFERENCES projects(id_project) ON DELETE CASCADE
);
--
CREATE TABLE publication_authors(
    id_project TEXT NOT NULL,
    id_user TEXT NOT NULL,
    author_position INT NOT NULL,  -- 0-indexed position of the user in the author list
    confirmed BIT DEFAULT 0,
    CONSTRAINT author_key PRIMARY KEY(id_project, id_user),
    FOREIGN KEY (id_project) REFERENCES projects(id_project),
    FOREIGN KEY (id_user) REFERENCES users_auth(id_user)
);

CREATE TABLE publication_external_authors(
	id_project TEXT NOT NULL,
	email TEXT NOT NULL,
	author_position INTEGER NOT NULL,
	confirmed BIT DEFAULT 0,
	FOREIGN KEY (id_project) REFERENCES projects(id_project),
	CONSTRAINT author_key PRIMARY KEY(id_project, email));

CREATE TABLE external_author_info(
	id_project TEXT NOT NULL,
	email TEXT NOT NULL,
	name TEXT NOT NULL,
	confirmed_datause BIT DEFAULT 0,
	confirmed_phi BIT DEFAULT 0,
	creation_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
	ip TEXT NOT NULL,
	FOREIGN KEY (id_project) REFERENCES projects(id_project),
	CONSTRAINT author_key PRIMARY KEY(id_project, email));

CREATE TABLE external_author_affiliation(
	id_project TEXT NOT NULL,
	email TEXT NOT NULL,
	affiliation TEXT NOT NULL);

CREATE TABLE project_invitations(
    id_invitation TEXT NOT NULL,
    id_invitee TEXT NOT NULL,  -- user being invited
    id_inviter TEXT NOT NULL,  -- user inviting
    id_project TEXT NOT NULL,  -- id of project to which the user is being invited
    time_sent TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (id_project) REFERENCES projects(id_project) ON DELETE CASCADE,
    CONSTRAINT invite_key PRIMARY KEY(id_invitee, id_project));

CREATE TABLE publications(
    id_publication VARCHAR(36) UNIQUE NOT NULL,
    id_project VARCHAR(36) NOT NULL,
    title VARCHAR(200) NOT NULL,
    description TEXT NOT NULL,
    creation_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    status VARCHAR(50) DEFAULT 'private',  -- if not specified, keep private
    hold BOOLEAN DEFAULT FALSE,
    id_user_contact VARCHAR(36),  -- id of the user to contact
    FOREIGN KEY (id_project) REFERENCES projects(id_project),
    FOREIGN KEY (id_user_contact) REFERENCES  users_auth(id_user));

CREATE TABLE project_papers(
    id_project VARCHAR(36) NOT NULL,
    doi VARCHAR(100) NOT NULL,
    url VARCHAR(200),
    title VARCHAR(200),
    year VARCHAR(5),
    author_str VARCHAR(200),
    num_authors INTEGER DEFAULT 0,
    FOREIGN KEY (id_project) REFERENCES projects(id_project),
    CONSTRAINT project_key PRIMARY KEY(id_project, doi)
);

CREATE TABLE publications_shortid(
    id_publication VARCHAR(36) UNIQUE NOT NULL,
    id_publication_short VARCHAR(5) UNIQUE NOT NULL,
    FOREIGN KEY (id_publication) REFERENCES publications(id_publication));

CREATE TABLE publications_analysis(
    id_publication VARCHAR(36) NOT NULL,
    id_analysis VARCHAR(36) NOT NULL,
    FOREIGN KEY (id_publication) REFERENCES publications(id_publication),
    FOREIGN KEY (id_analysis) REFERENCES analysis(id_analysis));

CREATE TABLE publications_data(
    id_publication VARCHAR(36) NOT NULL,
    id_data VARCHAR(36) NOT NULL,
    FOREIGN KEY (id_publication) REFERENCES publications(id_publication),
    FOREIGN KEY (id_data) REFERENCES data(id_data));

CREATE TABLE publications_results(
    id_publication VARCHAR(36) NOT NULL,
    id_result VARCHAR(36) NOT NULL,
    FOREIGN KEY (id_publication) REFERENCES publications(id_publication),
    FOREIGN KEY (id_result) REFERENCES results(id_result));

CREATE TABLE publications_authors(
    id_publication VARCHAR(36) NOT NULL,
    id_user TEXT NOT NULL,
    author_position INT NOT NULL,  -- 0-indexed position of the user in the author list
    confirmed BIT DEFAULT 0,
    CONSTRAINT author_key PRIMARY KEY(id_publication, id_user),
    FOREIGN KEY (id_publication) REFERENCES publications(id_publication),
    FOREIGN KEY (id_user) REFERENCES users_auth(id_user));

CREATE TABLE publications_authors_affiliation (
    id_publication VARCHAR(36) NOT NULL,
    id_user TEXT NOT NULL,
    affiliation TEXT NOT NULL,
    affiliation_order INTEGER DEFAULT 0,
    CONSTRAINT author_aff PRIMARY KEY (id_publication, id_user, affiliation),
    FOREIGN KEY (id_publication) REFERENCES publications(id_publication),
    FOREIGN KEY (id_user) REFERENCES users_auth(id_user));

CREATE TABLE publications_external_authors(
	id_publication VARCHAR(36) NOT NULL,
	email TEXT NOT NULL,
	author_position INTEGER NOT NULL,
	confirmed BIT DEFAULT 0,
	FOREIGN KEY (id_publication) REFERENCES publications(id_publication),
	CONSTRAINT author_key PRIMARY KEY(id_publication, email));

CREATE TABLE publications_external_authors_info(
	id_publication VARCHAR(36) NOT NULL,
	email TEXT,
	name TEXT,
    author_position INTEGER NOT NULL,
    confirmed_association BOOLEAN DEFAULT FALSE,
	confirmed_datause BOOLEAN DEFAULT FALSE,
	confirmed_phi BOOLEAN DEFAULT FALSE,
	creation_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
	FOREIGN KEY (id_publication) REFERENCES publications(id_publication),
	CONSTRAINT author_key PRIMARY KEY(id_publication, email));

CREATE TABLE publications_external_author_affiliation(
	id_publication VARCHAR(36) NOT NULL,
	email TEXT NOT NULL,
	affiliation TEXT NOT NULL,
    affiliation_order INTEGER DEFAULT 0);

CREATE TABLE publications_peer_review(
    id_publication VARCHAR(36) NOT NULL,
    peer_password TEXT NOT NULL,  -- password supplied to journal/reviewers
    FOREIGN KEY (id_publication) REFERENCES publications(id_publication) ON DELETE CASCADE
);

CREATE TABLE default_analysis(
    id_analysis VARCHAR(36) PRIMARY KEY NOT NULL,
    FOREIGN KEY (id_analysis) REFERENCES analysis(id_analysis)
)
