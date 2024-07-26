from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Flask, session, current_app, send_from_directory,
    jsonify, send_file
)


bp = Blueprint("tutorial", __name__, url_prefix="/tutorial")


@bp.route("/", methods=["GET"])
def landing():
    return render_template("tutorial/landing.html")

@bp.route("/pipeline_designer", methods=["GET"])
def pipeline_designer():
    return render_template("tutorial/pipeline_designer.html")