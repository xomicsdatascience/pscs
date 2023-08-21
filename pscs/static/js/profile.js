async function resetPassword(name_user){
    let msg = {"name_user": name_user};
    await fetch("/auth/reset/", {
        method: "POST",
        headers: {
            "Accept": "application/json",
            "Content-Type": "application/json"
        },
        body:JSON.stringify(msg)
    })
        .then(response => {return response.json()})
        .then(response => {console.log(response.url); window.location.href = response.url});
}

async function getUserAndReset(elId){
    let nameEl = document.getElementById(elId);
    let name_user = nameEl.value;
    await resetPassword(name_user);
}