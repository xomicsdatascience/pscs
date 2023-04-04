from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for, Flask, session, current_app, send_from_directory
)
from pscs.auth import admin_required
from pscs.db import get_db, get_unique_value_for_field
from markdown import markdown
bp = Blueprint("admin", __name__, url_prefix="/admin")


@bp.route("/post_creator", methods=["GET", "POST"])
@admin_required
def posts():
    """
    Manages the post creation page. These posts are stored in the database and rendered on the homepage.
    GET: Render page.
    POST: Expects form data with "title" and "content". Replaces HTML < with &lt; and stores the contents in the DB.
    """
    if request.method == "GET":
        return render_template("admin/new_posts.html")
    elif request.method == "POST":
        # New post is being submitted;
        form_data = request.form
        post_title = form_data["title"]
        post_content = form_data["content"]
        if "<script" in post_content:
            flash("Invalid input; post content can't contain HTML script tags")
            return redirect(url_for("admin.posts"))
        post_content.replace("<", "&lt;")

        # Put into database
        db = get_db()
        id_post = get_unique_value_for_field(db, "id_post", "posts")
        db.execute("INSERT INTO posts (id_post, title, body) VALUES (?, ?, ?)", (id_post, post_title, post_content))
        db.commit()
        flash("Post submitted")
        return redirect(url_for("pscs.index"))


@bp.route("/rendermd", methods=["POST"])
@admin_required
def rendermd():
    """
    Renders text into Markdown; expects JSON object keyed with "title" and "content"; only content is rendered.
    POST: Render the data in request.json["content"]. Warn the user if the post would get rejected on submission.
    """
    if request.method == "POST":
        # post_title = request.json["title"]
        post_content = request.json["content"]
        # before any markdown; this should not contain any "<"  (">" is used for block quotes in Markdown)
        post_content.replace("<", "&lt;")
        post_content = markdown(post_content)
        if "<script" in post_content:
            return {"title": request.json["title"], "content": "WARNING: Post can't contain &lt;script&gt; tags."}
        post = {"title": request.json["title"], "content": post_content}
        return post
