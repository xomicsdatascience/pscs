// JavaScript related to admins adding/editing posts on the newsfeed.
let idleTimer;
let idleTimeout = 1000;  // time in ms for when user is not typing

// document.addEventListener("keyup", waitForIdle);
document.addEventListener("keyup", waitForIdle);
function waitForIdle(){
    clearTimeout(idleTimer);
    idleTimer = setTimeout(renderMarkdown, idleTimeout);
}

async function renderMarkdown(){
    // sends the text in the content box to the server to have it render it and return
    let postTitleEl = document.getElementById("title");
    let postTitle = postTitleEl.value;

    let postContentEl = document.getElementById("content");
    let postContent = postContentEl.value;

    let postObj = {"title": postTitle, "content": postContent};

    let response = await fetch("/admin/rendermd", {
        method: "POST",
        headers: {
            "Accept": "application/json",
            "Content-Type": "application/json"
        },
        body: JSON.stringify(postObj)
    });
    let data = await response.json();

    // data obtained; set preview boxes
    let previewTitleEl = document.getElementById("previewTitle");
    let previewContentEl = document.getElementById("previewContent");
    previewTitleEl.innerHTML = data["title"];
    previewContentEl.innerHTML = data["content"];


}

