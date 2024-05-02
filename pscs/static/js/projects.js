/* ABOUT
This file contains JavaScript for the projects page of PSCS, mainly
the dynamic table creation based on analysis inputs.
*/
// consts
const INPUTDROPDOWNCLASS = "fileinputs"  // class name for SELECT element
const INPUTDROPDOWNPREFIX = "fileinput"  // id prefix for SELECT element
const IDSEP = '-'


lastTable = null;
function createTable(projectFiles, analysisInputs, buttonRunId='buttonRun'){
  // Create dropdown for project files
  drop = createDataDropdown(projectFiles);
  if(lastTable != null){
    lastEl = document.getElementById(lastTable);
    if(lastEl != null) {
      lastEl.remove();
    }
  }
  // Create row for each input node
  tbl = document.createElement('table');
  tbl.id = 'open';
  lastTable = tbl.id;
  for(var inp in analysisInputs){
    tr = tbl.insertRow();
    td = tr.insertCell();
    td.appendChild(document.createTextNode(analysisInputs[inp]));
    td = tr.insertCell();
    dropClone = drop.cloneNode(true);
    dropClone.class = INPUTDROPDOWNCLASS;
    dropClone.id = INPUTDROPDOWNPREFIX + IDSEP + inp.toString();
    td.appendChild(dropClone);
  }
  buttonRun = document.getElementById(buttonRunId);
  buttonRun.parentNode.insertBefore(tbl, buttonRun.nextSibling);
}


function createDataDropdown(files){
// Creates a SELECT menu based on the input "files" object; keys are used for the option.value, and values are used
// as the text.
  drop = document.createElement("select");
  drop.class = INPUTDROPDOWNCLASS;
  drop.id = INPUTDROPDOWNPREFIX + IDSEP + 'master';
  drop.style.maxWidth = '100px';
  for(var file in files){
    option = document.createElement("option");
    option.innerHTML = files[file];
    option.value = file;
    drop.appendChild(option);
  }
  return drop;
}

async function executePipeline(elId){
  let selectAnalysis = document.getElementById(elId);
  let idAnalysis = selectAnalysis.value;
  if(idAnalysis === "default"){
    return
  }
  // pause until we're done doing things
  pauseRun("buttonRun");

  let query = 'select[id^=' + INPUTDROPDOWNPREFIX + ']';
  let inputFiles = document.querySelectorAll(query);
  let filePaths = {};
  for(let i=0; i<inputFiles.length; i++){
    // .value has data id
    let filePath = inputFiles[i].value;
    if(filePath === ""){
      // no data associated with node
      unpauseRun("buttonRun");
      return
    }
    // .id contains the corresponding node's id
    let nodeId = parseFileInputId(inputFiles[i].id);
    filePaths[nodeId] = filePath;
  }
  let pipelineSummary = new Object();
  pipelineSummary['id_analysis'] = idAnalysis;
  pipelineSummary['file_paths'] = filePaths;

  await fetch("/run_analysis", {
    method: "POST",
    headers: {
      "Accept": "application/json",
      "Content-Type": "application/json"},
    body: JSON.stringify(pipelineSummary)
    }).then(response => {
      return response.json();
    })
      .then(data => {
        if(data.hasOwnProperty("submit_success") && data.submit_success === 0){
          alert("Job submission failed: " + data.submit_status);
        }
        else{
          alert("Job submitted!");
        }
      });
  unpauseRun("buttonRun")
}

let pauseStartCursor = null;
function pauseRun(elementIdToDisable){
  let disEl = document.getElementById(elementIdToDisable);
  pauseStartCursor = document.body.style.cursor;
  disEl.disabled = true;
  document.body.style.cursor = "wait";
}

function unpauseRun(elementIdToEnable){
  let disEl = document.getElementById(elementIdToEnable);
  disEl.disable = false;
  document.body.style.cursor = pauseStartCursor;
}

function startDeletion(id_data, name_data, file_hash){
  // Asks user to confirm their choice, then sends a request for the data to be deleted.
  const confirmationText = "Confirm deletion of:\n" + name_data + "\nwith hash\n" + file_hash;
  if(confirm(confirmationText)){
    let data_spec = new Object();
    data_spec['deleteData'] = id_data;
    fetch(window.location.href,{
      method: "POST",
      headers: {
          "Accept": "application/json",
          "Content-Type": "application/json"},
      body: JSON.stringify(data_spec)})
        .then(response => {window.location.href = response.url});
  }
}

function parseFileInputId(id){
  first_sep_idx = id.indexOf(IDSEP);
  return id.substr(first_sep_idx+1);
}

function delTable(id) {
  if(id == null){
    return;
  }
  tbl = document.getElementById(id);
  tbl.remove();
  lastTable = null;
  return
}

function validateRename(renameEl){
  if(renameEl.value.length > 0){
    return true
  }
  else{
    renameEl.setCustomValidity('New project name must be at least 1 character.');
    renameEl.reportValidity();
    return false
  }
}

function renameProject(){
  renameEl = document.getElementById('inputRename');
  if(validateRename(renameEl)){
    var rename_spec = new Object();
    rename_spec['newName'] = renameEl.value;
    fetch(window.location.href,{
    method: "POST",
    headers: {
        "Accept": "application/json",
        "Content-Type": "application/json"},
    body: JSON.stringify(rename_spec)
    }).then(response => {window.location.href = response.url});

  }
}

function deleteProject(projectName){
  confirmationText = "Delete project " + projectName + "?"
  if(confirm(confirmationText)){
    var delete_spec = new Object();
    delete_spec['delete'] = true;
    fetch(window.location.href, {
    method: "POST",
    headers: {
        "Accept": "application/json",
        "Content-Type": "application/json"},
    body: JSON.stringify(delete_spec)
    }).then(response => {window.location.href = response.url});
  }
}

function inviteUser(){
   let adduser_spec = {};
   let userEl = document.getElementById('inviteUser');
   adduser_spec['inviteUser'] = userEl.value;
   fetch(window.location.href, {
   method: "POST",
   headers: {
        "Accept": "application/json",
        "Content-Type": "application/json"},
   body: JSON.stringify(adduser_spec)
   }).then(response => {window.location.href = response.url});
}

function getTabInfo(tab){
  let idProject = window.location.href;
  idProject = idProject.split("/");
  idProject = idProject[idProject.length-1];
  fetch(idProject + '/tabs/' + tab)
      .then(response => response.text())
      .then(html => document.getElementById("tab_content").innerHTML = html)
      .catch(error => console.error("Error: ", error))
}

function displayResult(file_path){
  const container = document.getElementById("container_results");
  container.innerHTML = "<img src='/" + file_path + "'>" ;
  return
}

function viewLog(id_job){
  fetch(window.location.href + "/logs/" + id_job)
      .then(response => response.json())
      .then(data => setLogBoxes(data["stdout"], data["stderr"]))
  return
}

function setLogBoxes(stdout_data, stderr_data){
  const stdout_text = document.getElementById("stdout_text");
  const stderr_text = document.getElementById("stderr_text");
  stdout_text.value = stdout_data;
  stderr_text.value = stderr_data;
  const log_display = document.getElementById("log_display");
  log_display.style.display = "flex";
  return
}

function toggleDiv(divId){
    let div = document.getElementById(divId);
    if(div.style.display === "none"){
        div.style.display = "block";
    }
    else{
        div.style.display="none";
    }
}

function toggleClassDisplay(className){
  let elements = document.getElementsByClassName(className);
  for(let el of elements){
    if(el.style.display === "none"){
      el.style.display = "block";
    }
    else{
      el.style.display = "none";
    }
  }
}

function copyToClipboard(text){
  navigator.clipboard.writeText(text)
      .then(() => {})
      .catch(err => {console.error("Unable to copy text to clipboard: ", err)});
}



function notifyCopy(event){
  tempToolTip(event, "ID copied to clipboard.");
}

function tempToolTip(event, tempText){
  let tooltip = document.createElement("span");
  tooltip.innerHTML = tempText;
  tooltip.style.position = "fixed";
  tooltip.style.top = (event.clientY + 10) + "px";
  tooltip.style.left = (event.clientX - 20) + "px";
  tooltip.style.background = "#ddd";
  tooltip.style.border = "2px solid #000"
  tooltip.style.zIndex = 1000;
  tooltip.style.userSelect = "none";
  document.body.appendChild(tooltip);
  setTimeout(function(){tooltip.remove();}, 2000);
}

document.addEventListener('click', function(e) {
  let target = e.target;
  if(target.classList.contains("tab-button")){
    // get all tabs
    let tabs = document.getElementsByClassName("tab-button");
    for(let t of tabs) {
      t.classList.remove("tab-selected");
    }
    target.classList.add("tab-selected");
  }
}, false);


async function downloadPublishedResult(selectElID, warningElID){
  let sel = document.getElementById(selectElID);
  let analysis_obj = {"id_analysis": sel.options[sel.selectedIndex].value};
  let dest_url = window.location.pathname + "/request_results";
  await fetch(dest_url, {
    method: "POST",
    headers: {
      "Accept": "application/json",
      "Content-Type": "application/json"
    },
    body: JSON.stringify(analysis_obj)
  })
      .then(response => {
            let contentDisp = response.headers.get("Content-Disposition");
            let warningEl = document.getElementById(warningElID);
            if(response.status === 429){
                warningEl.style.display = "block";
                return
            }
            else{
                warningEl.style.display = "none";
            }
            let filename = contentDisp.split("filename=")[1];
            response.blob().then(blob => {
                let url = window.URL.createObjectURL(blob);
                let pretendLink = document.createElement("a");
                pretendLink.href = url;
                pretendLink.download = filename;
                pretendLink.click();
            });
        })
}

async function downloadResult(selectElID, warningElID){
  let sel = document.getElementById(selectElID);
  let analysis_obj = {"id_analysis": sel.options[sel.selectedIndex].value};
  let url = window.location.pathname;
  let urlcomp = url.split("/");
  let idx = urlcomp.indexOf("project");
  let project_url = urlcomp.slice(0, idx+2).join("/") + "/results_request"
  let response = await fetch(project_url, {
    method: "POST",
    headers: {
      "Accept": "application/json",
      "Content-Type": "application/json"
    },
    body: JSON.stringify(analysis_obj)
  })
  .then(response => {
            let contentDisp = response.headers.get("Content-Disposition");
            let warningEl = document.getElementById(warningElID);
            if(response.status === 429){
                warningEl.style.display = "block";
                return
            }
            else{
                warningEl.style.display = "none";
            }
            let filename = contentDisp.split("filename=")[1];
            response.blob().then(blob => {
                let url = window.URL.createObjectURL(blob);
                let pretendLink = document.createElement("a");
                pretendLink.href = url;
                pretendLink.download = filename;
                pretendLink.click();
            });
        })
  return
}
