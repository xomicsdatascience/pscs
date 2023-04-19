/* ABOUT
This file contains JavaScript for the projects page of PSCS, mainly
the dynamic table creation based on analysis inputs.
*/
// consts
const INPUTDROPDOWNCLASS = "fileinputs"  // class name for SELECT element
const INPUTDROPDOWNPREFIX = "fileinput"  // id prefix for SELECT element
const FIELDSEP = '-'


lastTable = null;
function createTable(projectFiles, analysisInputs){
  // Create dropdown for project files
  drop = createDataDropdown(projectFiles);
  if(lastTable != null){
    lastEl = document.getElementById(lastTable);
    lastEl.remove();
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
    dropClone.id = INPUTDROPDOWNPREFIX + FIELDSEP + inp.toString();
    td.appendChild(dropClone);
  }
  buttonRun = document.getElementById('buttonRun');
  buttonRun.parentNode.insertBefore(tbl, buttonRun.nextSibling);
}


function createDataDropdown(files){
// Creates a SELECT menu based on the input "files" object; keys are used for the option.value, and values are used
// as the text.
//  console.log(files);
  drop = document.createElement("select");
  drop.class = INPUTDROPDOWNCLASS;
  drop.id = INPUTDROPDOWNPREFIX + FIELDSEP + 'master';
  drop.style.maxWidth = '100px';
  for(var file in files){
    option = document.createElement("option");
    option.innerHTML = files[file];
    option.value = file;
    drop.appendChild(option);
  }
  return drop;
}

async function executePipeline(){
  // change cursor to wait
  let btn = document.getElementById("buttonRun");
  btn.disabled = true;
  let start_cursor = document.body.style.cursor
  document.body.style.cursor = "wait";
  query = 'select[id^=' + INPUTDROPDOWNPREFIX + ']';
  inputFiles = document.querySelectorAll(query);
  var filePaths = new Object();
  for(var i=0; i<inputFiles.length; i++){
    // .value has data id
    filePath = inputFiles[i].value;
    // .id contains the corresponding node's id
    console.log(inputFiles[i].id);
    nodeId = parseFileInputId(inputFiles[i].id);
    filePaths[nodeId] = filePath;
  }
  selectAnalysis = document.getElementById('analysis');
  idAnalysis = selectAnalysis.value;

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
      let data = response.json();
      return data;
    })
      .then(data => {
        console.log(data);
        if(data.hasOwnProperty("submit_success") && data.submit_success === 0){
          alert("Job submission failed: " + data.submit_status);
        }
        else{
          alert("Job submitted!");
        }
      });
  btn.disabled = false;
  document.body.style.cursor = start_cursor;
}

function startDeletion(id_data, name_data){
  // Asks user to confirm their choice, then sends a request for the data to be deleted.
  confirmationText = "Confirm deletion of:\n" + name_data
  if(confirm(confirmationText)){
    var data_spec = new Object();
    data_spec['deleteData'] = id_data;
    fetch(window.location.href,{
    method: "POST",
    headers: {
        "Accept": "application/json",
        "Content-Type": "application/json"},
    body: JSON.stringify(data_spec)
    }).then(response => {window.location.href = response.url});
  }
}

function parseFileInputId(id){
  first_sep_idx = id.indexOf(FIELDSEP);
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
  confirmationText = "Delete project " + projectName
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

function addUser(){
   var adduser_spec = new Object();
   userEl = document.getElementById('addUser');
   adduser_spec['addUser'] = userEl.value;
   fetch(window.location.href, {
   method: "POST",
   headers: {
        "Accept": "application/json",
        "Content-Type": "application/json"},
   body: JSON.stringify(adduser_spec)
   }).then(response => {window.location.href = response.url});
}

