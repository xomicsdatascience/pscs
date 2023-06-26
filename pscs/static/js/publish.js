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

function getCheckboxIds(checkboxClass){
    // Given a checkbx class, returns an array of the ID of the checkboxes of that class that are selected.
    let checkboxList = document.querySelectorAll("[class=" + checkboxClass + "]");
    const onlyChecked = Array.from(checkboxList).filter(checkboxList => checkboxList.checked);
    return onlyChecked.map(onlyChecked => onlyChecked.id);
}

function getAuthorInfo(authorSelectId="authorlist"){
    // Returns author objects, including external authors.
    let authlist = Array.from(document.getElementById(authorSelectId).options);
    // Extract only relevant info; id or name + email
    let authobj = [];
    for(let au of authlist){
        if(au.email == null){
            authobj.push({"id": au.id, "external": false});
        }
        else{
            authobj.push({"email": au.email, "external": true});
        }
    }
    return authobj
}

function addExternalAuthor(){
    // Adds the external author to the list; will get parsed later
    const authorSel = document.getElementById("authorlist");
    const emailEl = document.getElementById("external_email");
    let newOpt = document.createElement("option");
    if(emailEl.value === null || emailEl.value.length === 0){
        return;
    }
    newOpt.email = emailEl.value;
    newOpt.innerText = emailEl.value;
    authorSel.options.add(newOpt);
    emailEl.value = "";
    return;
}

function removeAuthor(selectId){
    // Removes the selected option from the <select> element with the input ID
    const selectEl = document.getElementById(selectId);
    if(selectEl.options.length === 1){
        selectEl.setCustomValidity("There must be at least one author.");
        selectEl.reportValidity();
        return;
    }
    let opt = selectEl.selectedOptions[0];
    if(!opt.hasOwnProperty("email")  || opt.email.length === 0){
        selectEl.setCustomValidity("PSCS users associated with the project can't be removed as authors.");
        selectEl.reportValidity();
        return;
    }
    selectEl.selectedOptions[0].remove();
}

async function parsePublicationPage(pubtype) {
    // Parses the publication page, creates a JSON object with the relevant info, then sends it to the server.
    let publicationSummary = {};
    publicationSummary["publication_type"] = pubtype;

    const confirmationEl = document.getElementById(pubtype);
    publicationSummary["confirmation"] = confirmationEl.value;

    // Get author ids, preserving order
    publicationSummary["authorlist"] = getAuthorInfo("authorlist");

    // Get analysis IDs
    publicationSummary["analyses"] = getCheckboxIds("analysisCheckbox");

    // Get data IDs
    publicationSummary["data"] = getCheckboxIds("dataCheckbox");

    await fetch(window.location.href, {
        method: "POST",
        headers: {
            "Accept": "application/json",
            "Content-Type": "application/json"
        },
        body: JSON.stringify(publicationSummary)
    })
        .then(response => {
            return response.json();
        })
        .then(data => {
            if (data.hasOwnProperty("publish_success") && data.publish_success === 1) {
                window.location.href = data.url;
            } else {
                alert("Project publication has failed; " + data.success_message);
            }
        });
}
