function moveUp(selectId){
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

function getAuthorIds(authorSelectId="authorlist"){
    // Returns the ids of the authors specified by the user
    let authlist = document.getElementById(authorSelectId);
    return Array.from(authlist.options).map(authlist => authlist.id);
}

function makeHiddenInput(name="", value=""){
    // Creates a hidden input and sets its name/value. Used for appending data to HTML forms.
    let hiddenInput = document.createElement("input");
    hiddenInput.type = "hidden";
    hiddenInput.name = name;
    hiddenInput.value = value;
    return hiddenInput;
}

function appendInputToForm(form, name, value){
    // Creates a hidden input with the name `name` and value `value` and appends it to `form`.
    let hiddenInput = makeHiddenInput(name, value);
    form.appendChild(hiddenInput);
}

async function parsePublicationPage(pubtype) {
    // Parses the publication page, creates a JSON object with the relevant info, then sends it to the server.
    let publicationSummary = {};
    publicationSummary["publication_type"] = pubtype;

    const confirmationEl = document.getElementById(pubtype);
    publicationSummary["confirmation"] = confirmationEl.value;

    // Get author ids, preserving order
    publicationSummary["authorlist"] = getAuthorIds("authorlist");

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
