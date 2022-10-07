DROP TABLE IF EXISTS users;
DROP TABLE IF EXISTS projects;
DROP TABLE IF EXISTS data;

CREATE TABLE users (
    id TEXT UNIQUE NOT NULL PRIMARY KEY,
    username TEXT UNIQUE NOT NULL,
    password TEXT NOT NULL,
    time_created TIMESTAMP DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE projects (
    id TEXT UNIQUE NOT NULL,
    name TEXT NOT NULL,
    description TEXT,
    userid TEXT NOT NULL,
    role TEXT NOT NULL,
    time_created TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    CONSTRAINT proj_key PRIMARY KEY (id,userid),
    FOREIGN KEY (userid) REFERENCES users(id)
);

CREATE TABLE data (
    id INTEGER PRIMARY KEY AUTOINCREMENT,  -- this should be updated to have uuid
    userid TEXT NOT NULL,
    file_path NOT NULL,
    sample_count INTEGER NOT NULL,
    gene_count INTEGER NOT NULL,
    file_hash TEXT NOT NULL,
    FOREIGN KEY (userid) REFERENCES users(id)
);
