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
        .then(response => {window.location.href = response.url});
}

async function getUserAndReset(elId){
    let nameEl = document.getElementById(elId);
    let name_user = nameEl.value;
    await resetPassword(name_user);
}

async function updateName(elId){
    // Updates the user's name
    let newName = document.getElementById(elId).value;

    await fetch("/auth/updateName", {
        method: "POST",
        headers: {
            "Accept": "application/json",
            "Content-Type": "application/json"
        },
        body: JSON.stringify({"newName": newName})
    })
        .then(response => {return response.json()})
        .then(response => {window.location.href = response.url});
}

function moveUp(selectId){
    // Moves up the selected option in the <select> element with the input ID
    const selectEl = document.getElementById(selectId);
    let selectedIdx = selectEl.selectedIndex;
    if(selectedIdx === 0 || selectedIdx === -1){
        return
    }
    let options = selectEl.options;
    const optionToAdd = options[selectedIdx];
    options.remove(selectedIdx);
    options.add(optionToAdd, selectedIdx-1);
    selectEl.options = options;
}


function moveDown(selectId){
    // Moves down the selected option in the <select> element with the input ID
    const selectEl = document.getElementById(selectId);
    let selectedIdx = selectEl.selectedIndex;
    let options = selectEl.options;
    if(selectedIdx === -1 || selectedIdx+1 === options.length){
        return
    }
    const optionToAdd = options[selectedIdx]
    options.remove(selectedIdx);
    options.add(optionToAdd, selectedIdx+1);
    selectEl.options = options;
}

function addAffiliation(sourceText, destinationSelect){
    const sourceEl = document.getElementById(sourceText);
    const newAff = sourceEl.value;
    const selectEl = document.getElementById(destinationSelect);
    const option = document.createElement("option");
    option.value = newAff;
    option.text = newAff;
    selectEl.add(option);
    sourceEl.value = "";
}

async function saveAffiliations(selectId){
    let affiliationsOpts = Array.from(document.getElementById(selectId).options);
    let affiliations = affiliationsOpts.map(option => option.value);
    await fetch("/auth/saveAffiliations", {
        method: "POST",
        headers: {
            "Accept": "application/json",
            "Content-Type": "application/json"
        },
        body: JSON.stringify({"affiliations": affiliations})
    })
        .then(response => {return response.json()})
        .then(response => {window.location.href = response.url});
}