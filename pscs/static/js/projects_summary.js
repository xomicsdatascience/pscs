async function respondToInvitation(id_invitation, action){
    let msg = {
        "action": action,
        "id_invitation": id_invitation
    };
    let response = await fetch("/project/manage_invitation", {
        method: "POST",
        headers: {
            "Accept": "application/json",
            "Content-Type": "application/json"
        },
        body: JSON.stringify(msg)
    });
    let data = await response.json();
    if(data.hasOwnProperty("url")){
        window.location.href = data.url;
    }
}

async function acceptInvitation(id_invitation){
    await respondToInvitation(id_invitation, "accept");
}

async function rescindInvitation(id_invitation){
    await respondToInvitation(id_invitation, "rescind");
}

async function rejectInvitation(id_invitation){
    await respondToInvitation(id_invitation, "reject");
}