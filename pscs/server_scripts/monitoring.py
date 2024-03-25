import os
import argparse
import sqlite3
from datetime import datetime, timedelta
from pathlib import Path
import shutil
from pscs.messaging.mail import send_email

sql_time_format = '%Y-%m-%d %H:%M:%S'

def get_users_between_dates(db: sqlite3.Connection,
                            start_date: str,
                            end_date: str) -> list:
    """
    Returns a list of users created between the specified times.
    Parameters
    ----------
    db : sqlite3.Connection
        Connection to the database
    start_date : str
        Date for when to start searching. Expects SQL-compatible format: YYYY-MM-DD HH:MM:SS. Use the Python datetime
        module and format string to generate and perform calculations with the time.
    end_date : str
        Date for when to stop searching. Expects SQL-compatible format: YYYY-MM-DD HH:MM:SS. Use the Python datetime
        module and format string to generate and perform calculations with the time.
    Returns
    -------
    list[dict]
        List of dicts keyed by 'id_user' and 'name_user' for each account created between the timepoints.
    """
    user_list = db.execute("SELECT id_user, name_user FROM users_auth "
                           "WHERE creation_time_user BETWEEN DATE(?) AND DATE(?)",
                           (start_date, end_date)).fetchall()
    return [dict(u) for u in user_list]


def get_projects_between_dates(db: sqlite3.Connection,
                               start_date: str,
                               end_date: str) -> list:
    """
    Returns a list of projects created between the specified times.
    Parameters
    ----------
    db : sqlite3.Connection
        Connection to the database
    start_date : str
        Date for when to start searching. Expects SQL-compatible format: YYYY-MM-DD HH:MM:SS. Use the Python datetime
        module and format string to generate and perform calculations with the time.
    end_date : str
        Date for when to stop searching. Expects SQL-compatible format: YYYY-MM-DD HH:MM:SS. Use the Python datetime
        module and format string to generate and perform calculations with the time.

    Returns
    -------
    list[dict]
        List of dicts keyed by 'id_project' and 'name_project' for each project created between the timepoints.
    """
    project_list = db.execute("SELECT id_project, name_project FROM projects "
                              "WHERE creation_time_project BETWEEN DATE(?) AND DATE(?)",
                              (start_date, end_date)).fetchall()
    return [dict(p) for p in project_list]


def get_published_projects_between_dates(db: sqlite3.Connection,
                                         start_date: str,
                                         end_date: str) -> list:
    """Similar to get_projects_between_dates, but for published projects."""
    project_list = db.execute("SELECT id_project, name_project FROM projects "
                              "WHERE is_published = 1 AND "
                              "date_published BETWEEN DATE(?) AND DATE(?)",
                              (start_date, end_date)).fetchall()
    return [dict(p) for p in project_list]


def get_submitted_jobs_between_dates(db: sqlite3.Connection,
                                     start_date: str,
                                     end_date: str) -> list:
    """
    Returns a list of jobs submitted between the specified times.
    Parameters
    ----------
    db : sqlite3.Connection
        Connection to the database.
    start_date : str
        Date for when to start searching. Expects SQL-compatible format: YYYY-MM-DD HH:MM:SS. Use the Python datetime
        module and format string to generate and perform calculations with the time.
    end_date : str
        Date for when to stop searching. Expects SQL-compatible format: YYYY-MM-DD HH:MM:SS. Use the Python datetime
        module and format string to generate and perform calculations with the time.

    Returns
    -------
    list[dict]
        List of dicts keyed by 'id_job', 'id_user', and 'id_project' that were submitted between the specified
        timepoints.
    """
    submitted_jobs = db.execute("SELECT id_job, id_user, id_project FROM submitted_jobs "
                                "WHERE date_submitted BETWEEN DATE(?) AND DATE(?)",
                                (start_date, end_date)).fetchall()
    return [dict(j) for j in submitted_jobs]


def get_new_users_in_last_days(db: sqlite3.Connection,
                               days: int = 7) -> list:
    """Gets the new users for the last X days."""
    start_date = (datetime.now() - timedelta(days=days)).strftime(sql_time_format)
    now_date = datetime.now().strftime(sql_time_format)
    return get_users_between_dates(db, start_date=start_date, end_date=now_date)


def get_new_projects_in_last_days(db: sqlite3.Connection,
                                  days: int = 7) -> list:
    """Gets the new projects for the last X days."""
    start_date = (datetime.now() - timedelta(days=days)).strftime(sql_time_format)
    now_date = datetime.now().strftime(sql_time_format)
    return get_projects_between_dates(db, start_date=start_date, end_date=now_date)


def get_new_published_projects_in_last_days(db: sqlite3.Connection,
                                            days: int = 7) -> list:
    """Gets newly-published projects for the last X days."""
    start_date = (datetime.now() - timedelta(days=days)).strftime(sql_time_format)
    now_date = datetime.now().strftime(sql_time_format)
    return get_published_projects_between_dates(db, start_date=start_date, end_date=now_date)


def get_submitted_jobs_in_last_days(db: sqlite3.Connection,
                                    days: int = 7) -> list:
    """Gets the submitted jobs for the last X days."""
    start_date = (datetime.now() - timedelta(days=days)).strftime(sql_time_format)
    now_date = datetime.now().strftime(sql_time_format)
    return get_submitted_jobs_between_dates(db, start_date=start_date, end_date=now_date)


def get_percent_disk_usage(path: Path) -> float:
    total, used, _ = shutil.disk_usage(path)
    return round(100*used / total, 2)

def get_disk_usage_gb(path: Path) -> (float, float):
    total, used, _ = shutil.disk_usage(path)
    return round(used / (1024 ** 3), 2), round(total / (1024 ** 3), 2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True, help="Path to database", type=str)
    parser.add_argument("--data_dir", required=False, help="Path to where data is stored.", type=str, default=None)
    parser.add_argument("--email", required=False, help="Email address to which to send the update. If not specified, the update is printed to the console.", type=str, default=None)
    args = vars(parser.parse_args())
    if args["data_dir"] is None:
        args["data_dir"] = os.path.dirname(args["db"])

    db = sqlite3.connect(args["db"], detect_types=sqlite3.PARSE_DECLTYPES)
    db.row_factory = sqlite3.Row
    db.execute('PRAGMA foreign_keys = ON;')

    # get updates in the last week
    new_users = len(get_new_users_in_last_days(db, days=7))
    new_projects = len(get_new_projects_in_last_days(db, days=7))
    new_publications = len(get_new_published_projects_in_last_days(db, days=7))
    new_submitted_jobs = len(get_submitted_jobs_in_last_days(db, days=7))
    data_disk_usage_pct = get_percent_disk_usage(args["data_dir"])
    data_disk_usage_gb_used, data_disk_usage_gb_total  = get_disk_usage_gb(args["data_dir"])

    db_disk_usage_pct = data_disk_usage_pct
    db_disk_usage_gb_used, db_disk_usage_gb_total = data_disk_usage_gb_used, data_disk_usage_gb_total
    if args["data_dir"] != os.path.dirname(args["db"]):
        db_disk_usage_pct = get_percent_disk_usage(args["db"])
        db_disk_usage_gb_used, db_disk_usage_gb_total = get_disk_usage_gb(args["db"])
    db.close()
    if data_disk_usage_pct > 75:
        data_disk_warning = "Data disk usage is high; consider increasing capacity soon."
    else:
        data_disk_warning = ""
    if db_disk_usage_pct > 75:
        db_disk_warning = "DB disk usage is high; consider increasing capacity soon."
    else:
        db_disk_warning = ""

    today_date = datetime.now().strftime("%Y/%m/%d")
    next_udpate = (datetime.now() + timedelta(days=7)).strftime("%Y/%m/%d")
    from pscs.templates.misc import monitoring_update_template
    update_text = monitoring_update_template.format(today_date=today_date,
                                                    next_update=next_udpate,
                                                    new_users=new_users,
                                                    new_projects=new_projects,
                                                    new_publications=new_publications,
                                                    new_submitted_jobs=new_submitted_jobs,
                                                    data_disk_usage_gb_used=data_disk_usage_gb_used,
                                                    data_disk_usage_gb_total=data_disk_usage_gb_total,
                                                    data_disk_usage_pct=data_disk_usage_pct,
                                                    data_disk_warning=data_disk_warning,
                                                    db_disk_usage_gb_used=db_disk_usage_gb_used,
                                                    db_disk_usage_gb_total=db_disk_usage_gb_total,
                                                    db_disk_usage_pct=db_disk_usage_pct,
                                                    db_disk_warning=db_disk_warning)
    if args["email"] is not None:
        send_email([args["email"]], subject="PSCS Weekly Update", body=update_text)
    else:
        print(update_text)