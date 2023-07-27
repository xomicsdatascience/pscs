
async function resetPassword(name_user){
    let msg = {"name_user": name_user};
    await fetch("/auth/reset/", {
        method: "POST",
        headers: {
            "Accept": "application/json",
            "Content-Type": "application/json"
        },
        body:JSON.stringify(msg)
    }).then(response => {window.location.href = response.url});
}