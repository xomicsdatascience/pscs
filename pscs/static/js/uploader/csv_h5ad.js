// addQuantListeners();
// addObsListeners();
// addVarListeners();
addListeners(quantListener, ["quantities", "quant_transpose"]);
addListeners(obsListener, ["obs_input", "obs_transpose"]);
addListeners(varListener, ["var_input", "var_transpose"]);

if(typeof window.csvFeatures === "undefined"){
    window.csvFeatures = "";
}
if(typeof window.csvIndex === "undefined") {
    window.csvIndex = "";
}

function addListeners(listener, elemIdArray){
    let elem;
    for(let elemId of elemIdArray){
        elem = document.getElementById(elemId);
        elem.addEventListener("change", listener);
    }
}

function quantListener(event){
    const quantEl = document.getElementById("quantities");
    let file = quantEl.files[0];
    if (file){
        const doTranspose = document.getElementById("quant_transpose").checked;
        displayText("quant_index", "Loading...");
        displayText("protein_names", "Loading...");
        if(!doTranspose) {
            processCSV(file, displayQuant);
        }
        else{
            processCSVTranspose(file, displayQuant);
        }
    }
}

function obsListener(event){
    const obsEl = document.getElementById("obs_input");
    let file = obsEl.files[0];
    if (file){
        const doTranspose = document.getElementById("obs_transpose").checked;
        displayText("obs_features","Loading...");
        displayText("obs_index","Loading...");
        if(!doTranspose){
            processCSV(file, displayObs);
        }
        else{
            processCSVTranspose(file, displayObs);
        }
    }
}

function varListener(event){
    const varEl = document.getElementById("var_input");
    let file = varEl.files[0];
    if(file){
        const doTranspose = document.getElementById("var_transpose").checked;
        displayText("var_features", "Loading...");
        displayText("var_index", "Loading...");
        if(!doTranspose){
            processCSV(file, displayVar);
        }
        else{
            processCSVTranspose(file, displayVar);
        }
    }
}
function displayText(element_id, text){
    const el = document.getElementById(element_id);
    el.textContent = text;
}

function displayQuant(results){
    if(results != null) {
        // data is fully loaded
        let keys = Object.keys(results.data[0]);
        let idx = keys[0];
        let proteinNames = keys[1];
        for (let s of keys.slice(2)) {
            proteinNames = proteinNames + "\n" + s;
        }
        displayText("quant_index", idx);
        displayText("protein_names", proteinNames);
    }
    else{
        displayText("quant_index", window.csvIndex);
        displayText("protein_names", window.csvFeatures.trim());
    }
}

function displayObs(results){
    if(results != null){
        // data is fully loaded (not transposed)
        let keys = Object.keys(results.data[0]);
        let idx = keys[0];
        let featureNames = keys[1];
        for(let s of keys.slice(2)){
            featureNames = featureNames + "\n" + s;
        }
        displayText("obs_index", idx);
        displayText("obs_features", featureNames);
    }
    else{
        displayText("obs_index", window.csvIndex);
        displayText("obs_features", window.csvFeatures.trim());
    }
}

function displayVar(results){
    if(results != null){
        // data is fully loaded (not transposed)
        let keys = Object.keys(results.data[0]);
        let idx = keys[0];
        let featureNames = keys[1];
        for(let s of keys.slice(2)){
            featureNames = featureNames + "\n" + s;
        }
        displayText("var_index", idx);
        displayText("var_features", featureNames);
    }
    else{
        displayText("var_index", window.csvIndex);
        displayText("var_features", window.csvFeatures.trim());
    }
}

function displayPreview(results, indexId, featureId){
    if(results != null){
        let keys = Object.keys(results.data[0]);
        let idx = keys[0];
        let featureNames = keys[1];
        for(let s of keys.slice(2)){
            featureNames = featureNames + "\n" + s;
        }
        displayText(indexId, idx);
        displayText(featureId, featureNames);
    }
    else{
        displayText(indexId, window.csvIndex);
        displayText(featureId, window.csvFeatures.trim());
    }
}

function processCSV(file, processFunc){
    Papa.parse(file, {
        preview: 1,
        header: true,
        complete: processFunc,
        worker: true
    });
}

function processCSVTranspose(file, processFunc){
    window.csvIndex = "";
    window.csvFeatures = "";
    Papa.parse(file, {
        header: true,
        worker: true,
        step: getFirstEl,
        complete: processFunc
    })
}
function processQuant(file){
    Papa.parse(file, {
        preview: 1,
        header: true,
        complete: displayQuant,
        worker: true,
    });
}


function getFirstEl(results){
    if(window.csvIndex === ""){
        window.csvIndex = Object.keys(results.data)[0];
        return
    }
    let data =results.data[window.csvIndex];
    if(window.csvFeatures === ""){
        window.csvFeatures = data;
    }
    else{
        window.csvFeatures = window.csvFeatures + "\n" + data;
    }
}