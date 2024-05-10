// This JS contains code for manipulating nodes, loading and saving pipelines

let lastSelectedElement = null;  // contains reference to last selected element; mostly used to delete it
let nodeIds = [];  // list of IDs for the nodes currently on the page
let lastOpenPanel = null;  // contains reference to the last opened panel; mostly used to close it
let movingNode = false;
// constants
// - id prefix of different elements
const SEP = "-";  // separating different fields in id string
const FIELDSEP = ".";  // separates data within a field (e.g. 5.2 for node 5, area 2)
const NODE = "pscs";
const NODECLASS="pscsNode";
const CONNECTOR = "connector";
const IMAP = "imagemap";
const IAREA = "imagearea";
const IAREA_IN = IAREA + SEP + "in";
const IAREA_OUT = IAREA + SEP + "out";
const LABEL = "label";
const LABELCLASS = "pscsLabel"
const PARAMPANEL = "pscsParamPanel";
const INPUT = "input";
const PARAMNAME = "paramName";
const NODENAME = "pscsNodeName";
const EXPORTPANEL = "exportPanel";

// - constant values
const NODE_WIDTH = 60;  // width of image representing a node
const NODE_HEIGHT = 60;  // height of image representing a node
const XPUT_WIDTH = 20;  // size of input/output connections
const IAREA_MARGIN = 10;  // number of pixels between adjacent input/outputs
const DELKEY = 46;
const MACDELKEY=8;

const VERTICAL_POSITION_INITIAL = 200;  // initial vert position for new nodes (px)
const HORIZONTAL_POSITION_INITIAL = 200;  // initial horizontal position for new nodes (px)
const VERTICAL_MARGIN = 30;  // margin between new element and previous ones (px)
const HORIZONTAL_MARGIN = 30;  // margin between new element and previous ones (px)
const VERTICAL_SPACING = 100;
const HORIZONTAL_SPACING = 100;
const VERTICAL_JIGGLE = 10;
const HORIZONTAL_JIGGLE = 15;

let PADDING_LEVEL = 10;  // sidebar padding
let START_FONTSIZE = 20;  // sidebar fontsize for top-level module
let DECREASE_STEP_FONTSIZE = 4;  // decrease fontsize by this amount for every nested module
let MIN_FONTSIZE = 12;  // minimum font size for nested module
// const ISTR = 'param[thing]';
const PARAM_KEYWORD = 'param';
const PARAM_REGEX = new RegExp(`${PARAM_KEYWORD}\\[([^\\]]+)\\]`);
const PARAM_LISTSTRINGS = ["Collection", "Sequence", "list"];
// const result = ISTR.match(regex);

let nodeInfo = null;
getNodeInfo();

// document listeners
document.addEventListener("keydown", kbshortcut);
function kbshortcut(event){
    if(event.keyCode === DELKEY || event.keyCode === MACDELKEY) {
        // in case user is deleting something in a textbox:
        if(!(document.activeElement instanceof HTMLInputElement)) {
            if (lastSelectedElement != null) {
                lastSelectedElement.del();
                lastSelectedElement = null;
            }
            if (lastOpenPanel !== null){
                lastOpenPanel.remove();
                lastOpenPanel = null;
            }
        }
    }
}

function formatId(idPrefix, ...args){
// Formats an id string with the structure idPrefix + SEP + args[0] + SEP + [...]
  let retString = idPrefix;
    for(let idx=0; idx<args.length; idx++){
      retString += SEP + args[idx].toString();
    }
  return retString;
}

function extractIdNums(idStr){
// Extracts the number from the string idStr, expecting it to be of form idPrefix + SEP + num + SEP + [...].
// Returns: array, with each entry corresponding to the parts of the input IdStr.
  let idSplit = idStr.split(SEP).slice(1);  // remove id prefix
  // some entries will have FIELDSEP in them; return array at that pos
  let idArr = [];
  idSplit.forEach(element => {
      idArr.push(element.split(FIELDSEP));
  });
  return idArr;
}

function getUniqueId(idPrefix){
// Gets a unique ID of the form "idPrefix" + SEP + [int]
  // gather all elements that have the same prefix
  const prefixElements = document.querySelectorAll("[id^=" + idPrefix);
  let newNum = 0;
  // get ids of similar elements, identify likely number to create unique name
  if(prefixElements.length > 0){
    let prefixElementIds = [];
    prefixElements.forEach(element => prefixElementIds.push(element.id));
    prefixElementIds.sort();
    let lastNum = extractIdNums(prefixElementIds.slice(-1)[0])[0];
    newNum = parseInt(lastNum) + 1;
  }
  // in case something went weird with the sorting, increment the number until you get a number
  let tentativeId = formatId(idPrefix, newNum);
  let extantElement = document.getElementById(tentativeId);
  while(extantElement != null){
    newNum += 1;
    tentativeId = formatId(idPrefix, newNum);
    extantElement = document.getElementById(tentativeId);
  }
  return tentativeId;
}

function moveElement(moveX, moveY, transitionTime=0){
// moves an element on the page by moveX, moveY in the X- and Y- directions. Transition time determines how long the
// move should take.
// moveX: float - amount to move in the X direction
// moveY: float - amount to move in the Y direction
// transitionTime: float - amount of time in seconds that the transition should last
  // store original value to reset later
  const origTransititionDuration = this.style.transitionDuration;
  this.style.transitionDuration = transitionTime.toString() + "s";
  // update positions
  this.style.left = addToPx(this.style.left, moveX);
  this.style.top = addToPx(this.style.top, moveY);
  // reset transition duration once we're done moving
  this.addEventListener("transitionend", function resetTransitionDuration(){
    this.style.transitionDuration = origTransititionDuration;  //reset
    this.removeEventListener("transitionend", resetTransitionDuration);
  });
}

function moveElementToPos(posX, posY, transitionTime=0){
// moves an element to a specific position in px.
// posX: float - X position to move to
// posY: float - Y position to move to
// transitionTime: float - amount of time in seconds the move should take
  const origTransititionDuration = this.style.transitionDuration;
  this.style.transitionDuration = transitionTime.toString() + 's';
  this.style.left = posX.toString() + "px";
  this.style.top = posY.toString() + "px";
  this.addEventListener("transitionend", function resetTransitionDuration(){
    this.style.transitionDuration = origTransititionDuration;  //reset
    this.removeEventListener("transitionend", resetTransitionDuration);
  });
}

function addToPx(pxString, value){
  const pxValue = getValueFromPx(pxString);
  return (pxValue + value).toString() + "px";
}

function getValueFromPx(pxString){
  return parseFloat(pxString.split("px")[0]);
}

function createPscsNode(idNum, img, nodeData){
  // idNum: int - optional; if defined, will attempt to make a node with this id: NODE-idNum
  // processName: string - name of the function; used for initial labelText
  // module: string - name of the Python module that the process came from
  // params: string - JSON string containing process parameters of form: paramName: [paramType, paramValue]
  // pscsType: string - type of node: {"input", "output", "simo", "mimo"}
  // img: string - path to the image to use for this node
    let processName = nodeData["name"];
    let module = nodeData["module"];
    let params = nodeData["parameters"];
    let pscsType = nodeData["type"];
    // Get Id for node
    let pageEl = document.createElement("img");
    pageEl.class = NODE;
    pageEl.classList.add(NODE);
    pageEl.classList.add(NODECLASS);
    applyClassCSS(pageEl);
    let nodeId = null;
    // check if input ID is valid
    if(idNum != null && idNum.length !== 0){
      let tentativeId = formatId(NODE, idNum);
      let extantNodeWithId = document.getElementById(tentativeId);
      if(extantNodeWithId != null){
        return;
      }
      else{nodeId = tentativeId}
    }
    else{nodeId = getUniqueId(NODE)}
    pageEl.id = nodeId;
    pageEl.nodeId = nodeId;  // to ensure that this is exported
    // Assign name, module
    pageEl.procName = processName;
    pageEl.module = module;

    pageEl._requirements = structuredClone(nodeData["requirements"]);
    pageEl._effects = structuredClone(nodeData["effects"]);

    pageEl.requirements = structuredClone(pageEl._requirements);
    pageEl.effects = structuredClone(pageEl._effects);

    nodeIds.push(pageEl.id);  // add id to the list -- should probably remove this and use querySelectorAll instead
    pageEl.params = [];
    pageEl.important_parameters = [];
    pageEl.paramsValues = {};
    pageEl.paramsTypes = {};
    // store param names, values into pageEl
    for(let p of params){
        pageEl.params.push(p["name"]);
        pageEl.paramsTypes[p["name"]] = p["type"];
        pageEl.paramsValues[p["name"]] = p["default"];
    }

    if(nodeData !== null) {
        pageEl.important_parameters = nodeData["important_parameters"];
        pageEl.required_parameters = nodeData["required_parameters"];
        pageEl.num_inputs = nodeData["num_inputs"];
        pageEl.num_outputs = nodeData["num_outputs"];
    }
    pageEl.title = pageEl.params.toString();
    pageEl.pscsType = pscsType;
    if(!img.startsWith("/")){
        img = "/" + img;
    }
    pageEl.img = img; // for later loading
    pageEl.src = img;
    pageEl.srcConnectors = [];
    pageEl.dstConnectors = [];
    pageEl.coMove = [];  // list of things to move at the same time

    // Create label
    pageEl.labelText = pageEl.procName;  // user can set label name, but not procName

    // Create image map
    let [numInput, numOutput] = getNumInputOutputs(pageEl);
    if(pscsType === "input"){numInput = 0}
    else if(pscsType === "output"){numOutput=0}
    pageEl.imagemap = createImageMap(pageEl.id, numInput, numOutput);
    pageEl.useMap = "#" + pageEl.imagemap.id;
    // add element to page
    pageEl.style.zIndex = "3";
    pageEl.style.left = "0px";
    pageEl.style.top = "0px";
    addElementToContainer(pageEl);

    // attach methods
    pageEl.move = moveElement;
    pageEl.moveAll = moveAll;
    pageEl.onmousedown = grab;
    pageEl.label = createLabel(pageEl, pageEl.labelText, 0, NODE_HEIGHT);
    pageEl.del = del;
    pageEl.ondblclick = openPanel;
    pageEl.setInvalid = setInvalid;
    pageEl.resolveInteractions = resolveInteractions;

    // define methods
    let startX = 0;
    let startY = 0;
    function grab(event){
    // signals that the element has been clicked on and should follow the cursor
      event.preventDefault();
      event.target.style.cursor = "grabbing";
      if(lastSelectedElement !== null) {
          lastSelectedElement.coMove = lastSelectedElement.coMove.filter(el => el.id !== "selectionMarker");
          delMarker();
      }
      lastSelectedElement = pageEl;  // for deletion
        let selectionMarker = makeSelectionMarker();
        alignElementCenterToOther(selectionMarker, pageEl, [-15,0]);
        lastSelectedElement.coMove.push(selectionMarker);

      startX = event.clientX;  // where mouse started
      startY = event.clientY; // where mouse started
      document.onmouseup = release;  // for when user releases mouse
      document.onmousemove = drag;  // to move with cursor
      movingNode = true;  // signal that node is being moved; prevents dot from appearing
    }
    function release(event){
        // signals that the element should stop following the cursor
        event.preventDefault();
        lastSelectedElement.style.cursor = "grab";
        document.onmouseup = null;
        document.onmousemove = null;
        movingNode = false;  // signal that node is no longer being moved
        // expand container / prevent shrinking
        updateContainerSize();
        // limit node position to container
        let containerOffset = getContainerOffset();
        let pos = pageEl.getBoundingClientRect();
        let moveX = 0;
        let moveY = 0;
        if(pos.left < containerOffset[0]){
            moveX = containerOffset[0] - pos.left;
        }
        if(pos.top < containerOffset[1]){
            moveY = containerOffset[1] - pos.top;
        }
        if(moveX != 0 || moveY != 0){
            pageEl.moveAll(moveX, moveY);
        }
    }
    function drag(event){
    // signals that the element should be moved to follow cursor
      event.preventDefault();
      // get amount that the cursor has moved
      let mouseX = event.clientX;
      let mouseY = event.clientY;
      let moveX = mouseX - startX;
      let moveY = mouseY - startY;
      // set start position to mouse position
      startX = mouseX;
      startY = mouseY;
      // do move
      pageEl.moveAll(moveX, moveY);
    }
    function moveAll(moveX, moveY){
    // moves the element and all associated elements
      pageEl.move(moveX, moveY);
      pageEl.label.move(moveX, moveY);
      for(let idx=0; idx<pageEl.srcConnectors.length; idx++){
          let connector = document.getElementById(pageEl.srcConnectors[idx]);
          connector.moveSrc(moveX, moveY);
      }
      for(let idx=0; idx<pageEl.dstConnectors.length; idx++){
          let connector = document.getElementById(pageEl.dstConnectors[idx]);
          connector.moveDst(moveX, moveY);
      }
      for(let idx=0; idx<pageEl.coMove.length; idx++){
          pageEl.coMove[idx].move(moveX, moveY);
      }
    }
    function del(){
        // delete this node and all attached elements
        pageEl.label.del();
        for(let idx=pageEl.srcConnectors.length-1; idx>=0; idx--){
        // for(let idx=0; idx<pageEl.srcConnectors.length; idx++){
            // removing first element moves the array; keep deleting the first element until none are left
            document.getElementById(pageEl.srcConnectors[idx]).del();
        }
        for(let idx=pageEl.dstConnectors.length-1; idx>=0; idx--){
        // for(let idx=0; idx<pageEl.dstConnectors.length; idx++){
            document.getElementById(pageEl.dstConnectors[idx]).del();
        }
        for(let idx=pageEl.coMove.length-1; idx>=0; idx--){
            pageEl.coMove[idx].remove();
        }
        nodeIds = nodeIds.filter(id => id !== pageEl.nodeId);
        lastSelectedElement = null;
        let delMarker = document.getElementById("selectionMarker");
        if(delMarker !== null){
            delMarker.remove();
        }
        pageEl.remove();
    }
    function openPanel(){
        // Close other panel if it's already opened
        if (lastOpenPanel != null){
            lastOpenPanel.remove();
        }
        // Don't open a panel if node has no parameters
        if (pageEl.params.length === 0){
            return;
        }
        lastSelectedElement.coMove = lastSelectedElement.coMove.filter(el => el.id !== "selectionMarker");
        lastSelectedElement = null;  // prevent character deletion from user hitting 'delete'
        delMarker();
        const panel = document.createElement("div");
        panel.class = PARAMPANEL;
        panel.id = formatId(PARAMPANEL, extractIdNums(pageEl.id)[0]);
        panel.nodeId = pageEl.id;
        panel.describedNode = pageEl;
        applyClassCSS(panel);
        // Create table for user input
        let tbl = document.createElement('table');
        tbl.classList.add("pscsParamTable");
        // Name for display / input selection later
        let tr = tbl.insertRow();
        let td = tr.insertCell();
        td.appendChild(document.createTextNode('Node label'));
        td = tr.insertCell();
        let inp = document.createElement('input');
        inp.defaultValue = pageEl.labelText;
        inp.id = formatId(NODENAME, extractIdNums(pageEl.id)[0]);
        td.appendChild(inp);
        let displayParams;
        let displayValues = {};
        let displayTypes = {};
        let important_only = false;
        // If important parameters are listed, include only those
        if (pageEl.important_parameters !== null) {
            important_only = true;
            displayParams = pageEl.important_parameters;
            for (let idx in displayParams){
                let dName = displayParams[idx];
                displayValues[dName] = pageEl.paramsValues[dName];
                displayTypes[dName] = pageEl.paramsTypes[dName];
            }
        }
        else{
            important_only = false;
            displayValues = pageEl.paramsValues;
            displayTypes = pageEl.paramsTypes;
        }

        populatePanel(panel, displayValues, displayTypes, important_only);
        addElementToContainer(panel);
        lastOpenPanel = panel;
    }

    function populatePanel(panel, paramsValues, paramsTypes, important_only=false){
        const currentTable = panel.querySelector("table");
        if (currentTable !== null) {
            currentTable.remove();
        }
        const tbl = document.createElement("table");
        tbl.classList.add("pscsParamTable");
        let tr, td;
        let parVal, parType;
        for (const parName in paramsValues){
            parVal = paramsValues[parName];
            parType = paramsTypes[parName];

            tr = tbl.insertRow();
            td = tr.insertCell();
            td.appendChild(document.createTextNode(parName));

            td = tr.insertCell();
            const inp = document.createElement('INPUT');
            if (parVal != null){
                inp.defaultValue = parVal;
            }
            inp.id = formatId(PARAMNAME, parName);
            inp.parType = parType;
            td.appendChild(inp);
        }

        // Create save/close buttons
        tr = tbl.insertRow();
        td = tr.insertCell();
        let btn = document.createElement("button");
        btn.onclick = function(){closePanel(btn.panelId)};
        btn.panelId = panel.id;
        btn.innerHTML = "&times; Cancel"
        td.appendChild(btn);

        td = tr.insertCell();
        btn = document.createElement("button");
        btn.panelId = panel.id;
        const pageEl = panel.describedNode;
        if (important_only){
            // Add button that displays all parameters
            btn.onclick = function(){populatePanel(panel, pageEl.paramsValues, pageEl.paramsTypes, important_only=false)}
            btn.innerHTML = "Expand";
        }
        else{
            btn.onclick = function(){openPanel();}
            btn.innerHTML = "Shrink";
        }
        td.appendChild(btn);
        td = tr.insertCell();
        btn = document.createElement("button");
        btn.onclick = saveParams;
        btn.nodeId = pageEl.id;
        btn.panelId = panel.id;
        btn.innerHTML = "&#10003; Save"
        td.appendChild(btn);
        panel.appendChild(tbl);
    }

    function closePanel(panelId){
        let panel = document.getElementById(panelId);
        panel.remove();
    }
    function saveParams(event){
        // Save parameters that are in the panel and store them in the node.
        let btn = event.target;
        let panel = document.getElementById(btn.panelId);
        let tbl = panel.children[0];
        let params = {};
        for (let i=0, row; row = tbl.rows[i]; i++){
            for (let j=0, cell; cell = row.cells[j]; j++){
                if (cell.children.length > 0 && cell.children[0].id.includes(PARAMNAME)){  // only select boxes of this class
                    let txtBox = cell.children[0];
                    let param = txtBox.id.substring(PARAMNAME.length + 1, txtBox.id.length);
                    if (txtBox.value.length > 0){
                        params[param] = txtBox.value;
                    }
                    else if (txtBox.value.length === 0){
                        params[param] = null;
                    }
                }
                else if(cell.children.length > 0 && cell.children[0].id.includes(NODENAME)){
                    // display name for the node
                    let txtBox = cell.children[0];
                    let name = txtBox.value;
                    pageEl.updateLabel(name);
                }
            }
        }
        for(let p in params){
            pageEl.paramsValues[p] = params[p];
        }
        closePanel(panel.id);
    }
    pageEl.updateLabel = updateLabel;
    function updateLabel(text){
        if(text !== null) {
            pageEl.labelText = text;
        }
        pageEl.label.update(pageEl.labelText);
    }

    function setInvalid(text){
        // If the node is determined to be invalid, mark it with a dot that can be hovered over to identify the problem
        let dot = document.createElement("span");
        dot.classList.add("invalidMarker");
        dot.move = moveElement;
        dot.title = text;
        dot.style.opacity = "0%";
        addElementToContainer(dot);

        // Make centers align
        let dotWidth = 15;
        let dotHeight = 15;
        let parentRec = pageEl.getBoundingClientRect();
        let centerX = parentRec.left + parentRec.width/2;
        let centerY = parentRec.top + parentRec.height/2;
        const [contX, contY] = getContainerOffset();
        dot.style.left = (centerX - dotWidth/2-contX) + "px";
        dot.style.top = (centerY - dotHeight/2-contY) + "px";
        pageEl.coMove.push(dot);
        dot.style.transitionDuration = "1s";
        dot.style.transitionProperty = "opacity";
        dot.style.opacity = "100%";
    }

    function resolveInteractions(){
        // Resolves the parameterized interaction strings into strings with the node's values.
        this.requirements = structuredClone(this._requirements);
        this.effects = structuredClone(this._effects);
        for(let metaInteraction of [this.requirements, this.effects]) {
            for (let interaction of metaInteraction) {
                for (let key in interaction) {
                    if (interaction[key].length > 0) {
                        // let req_replace = req[key];
                        let interactionReplace = [];
                        let resolvedString = [];
                        for (let interactionValue of interaction[key]) {
                            resolvedString = resolveString(interactionValue, this.paramsValues, this.paramsTypes);
                            interactionReplace = interactionReplace.concat(resolvedString);
                        }
                        interaction[key] = interactionReplace;
                    }
                }
            }
        }
        let reqRemove = [];
        if(this.requirements.length > 1) {
            for (let idx0 in this.requirements) {
                let req0 = this.requirements[idx0];
                for (let req1 of this.requirements.slice(idx0+1)) {
                    let count0 = getInteractionCount(req0);
                    let count1 = getInteractionCount(req1);
                    let overlap = getInteractionOverlapCount(req0, req1);
                    if(count0 === count1 && count0 === overlap){
                        reqRemove.push(idx0);
                    }
                }
            }
        }
        for(let idx of reqRemove.reverse()){
            this.requirements.splice(idx, 1);
        }
    }

    function resolveString(s, paramsValue, paramsTypes){
        // Check for matching param string
        let formattedStringArray = formatInteractionString(s, paramsValue, paramsTypes);
        formattedStringArray = formattedStringArray.filter(el=>el!==null);
        return formattedStringArray
    }

    pageEl.cumulativeEffect = cumulativeEffect;
    function cumulativeEffect(firstCall=true){
        // determine the effect up to this node
        this.resolveInteractions();
        let cumulEffect = [];
        if(!firstCall){
            cumulEffect = structuredClone(this.effects);
        }
        let toCombine = [];
        for(let prev of getPrevNodes(this)){
            toCombine = toCombine.concat(prev.cumulativeEffect(false));
        }
        let toReturn = [];
        if(cumulEffect.length === 0){
            toReturn = structuredClone(toCombine);
        }
        else if(toCombine.length === 0){
            toReturn = structuredClone(cumulEffect);
        }
        else {
            for (let cumul of cumulEffect) {
                for (let comb of toCombine) {
                    toReturn.push(combineInteraction(cumul, comb));
                }
            }
        }
        if(toReturn.length > 0) {
            this.cumul = toReturn;
        }
        else{
            this.cumul = [{"obs": [], "obsm": [], "obsp": [], "var": [], "varm": [], "var_names": [], "uns": [], "layers": []}]
        }
        return toReturn;
    }
    pageEl.evaluateReady = evaluateReady;
    function evaluateReady(){
        let unmetRequirements = [];
        // Evaluates whether this node's requirements have been met by previous nodes.
        for(let requirementSet of this.requirements){  // check each set in requirements
            // Check whether all contents of this requirement is in one of the cumulative effects
            let unmetReqSet = {};
            let numUnmet = 0;
            for(let key in requirementSet){
                unmetReqSet[key] = [];
                for(let req of requirementSet[key]){
                    let isMet = false;
                    for(let eff of this.cumul){
                        if(eff[key].includes(req)){
                            isMet = true;
                            break
                        }
                    }
                    if(!isMet){
                        unmetReqSet[key].push(req);
                        numUnmet += 1;
                    }
                }
            }
            if(numUnmet > 0) {
                unmetRequirements.push(unmetReqSet);
            }
        }
        this.unmetRequirements = unmetRequirements;
        return this.unmetRequirements.length <= 0;
    }

    function combineInteraction(inter0, inter1){
        let newInter = {};
        for(let key in inter0){
            newInter[key] = inter0[key].concat(inter1[key]);
        }
        return newInter;
    }

    return pageEl;
}

function formatInteractionString(string, replacements, replacementTypes) {
    // Takes a string that contains keyed values determined by PARAM_REGEX
    let replacementArray = [];
    //go through keys that are in the string
    let parKey = string.match(PARAM_REGEX);
    while (parKey !== null) {
        let val = replacements[parKey[1]];
        if (val === null || val.length === 0) {
            return [null]
        }
        if(couldBeList(replacementTypes[parKey[1]])){
            // split values into list
            val = String(val).split(",").map(el=>el.trim());
            for(let v of val){
                let s = string.replace(parKey[0], v);
                replacementArray = replacementArray.concat(formatInteractionString(s, replacements, replacementTypes));
            }
            return replacementArray;
        }
        string = string.replace(parKey[0], val)
        parKey = string.match(PARAM_REGEX);
    }
    return [string];
}

function couldBeList(interactionString){
        for(let p of PARAM_LISTSTRINGS){
            if(interactionString.includes(p)){
                return true;
            }
        }
        return false;
}

function selectPackage(packageName, packageList){
    for(let p of packageList){
        if(p["modules"]["name"] === packageName){
            return p
        }
    }
    return null
}


function getNodeDataFromModule(nodeModule, nodeName, moduleInfo){
    let split = nodeModule.split(".");
    let sub = split.slice(1).join(".");
    let next = split[1];
    for (let m of moduleInfo["modules"]) {
        if (m["name"] === next) {
            if (sub.includes(".")) {
                return getNodeDataFromModule(sub, nodeName, m);
            } else {
                // bottom level; go through nodes
                for (let n of m["nodes"]) {
                    if (n["name"] === nodeName) {
                        return n
                    }
                }
            }
        }
    }
    return null
}

function updateContainerSize(){
    let container = getContainer();
    // get maximum position of nodes
    let nodes = document.querySelectorAll("[id^=" + NODE + SEP + "]");
    let maxRight = 0;
    let maxBottom = 0;
    nodes.forEach(node => {
        let rect = node.getBoundingClientRect();
        maxRight = Math.max(rect.right + NODE_WIDTH/4, maxRight);
        maxBottom = Math.max(rect.bottom + NODE_HEIGHT/4, maxBottom);
    });
    let containerOffset = getContainerOffset();
    container.style.minWidth = (maxRight - containerOffset[0]).toString() + "px";
    container.style.minHeight = (maxBottom - containerOffset[1]).toString() + "px";
}

function getNumInputOutputs(nodeEl){
    // gets the number of input/output ports
    return [nodeEl.num_inputs, nodeEl.num_outputs]
}

function createLabel(nodeEl, labelText, xPos=0, yPos=0){
  // creates a label and adds it to the node container at the specified position
  // nodeId: string - ID of the node this label is attacked to
  // labelText: string - text to use for label
    // create HTML element
    let labelEl = document.createElement("div");
    // create id
    const nodeNum = extractIdNums(nodeEl.id)[0];
    labelEl.id = formatId(LABEL, nodeNum);
    labelEl.class = LABELCLASS;
    applyClassCSS(labelEl);
    labelEl.pscsNode = nodeEl.id;
    labelEl.innerHTML = labelText;
    labelEl.style.left = xPos.toString() + "px";
    labelEl.style.top = yPos.toString() + "px";
    labelEl.style.position = "absolute";
    labelEl.onmousedown = nodeEl.onmousedown;
    labelEl.style.zIndex = nodeEl.style.zIndex;
    // attach methods
    labelEl.move = moveElement;
    labelEl.del = del;
    labelEl.update = update;
    addElementToContainer(labelEl);
    let beingRemoved = false;
    return labelEl;

    function del(){
        // Deleting the node causes the label to be deleted, causing a stack overflow.
        // This lets us get around it.
        if(beingRemoved){
            labelEl.remove();
        }
        else{
            beingRemoved = true;
            let attachedNode = getPscsNodeFromNum(nodeNum);
            attachedNode.del();
        }

    }
    function update(text){
        labelEl.innerHTML = text.toString();
    }
}

function createImageMap(nodeId, numInput=1, numOutput=1){
// creates an image map associated with nodeId.
// nodeId: string - id of the node with which to associate the imagemap
// numInput: int - number of inputs
// numOutput: int - number of outputs
  // get nodeNum
  let nodeNum = extractIdNums(nodeId)[0];
  let mapEl = document.createElement("map");
  mapEl.id = formatId(IMAP, nodeNum);
  mapEl.style.cursor = "grab";

  // create areas
  // inputs
  let input_linspace = linspace(0, NODE_HEIGHT, numInput+2);  // +2 and we'll drop the ends; this centers the ports
  input_linspace = input_linspace.slice(1,-1);
  createAreasForImageMap(mapEl, input_linspace,"in");
  let output_linspace = linspace(0, NODE_HEIGHT, numOutput+2);
  output_linspace = output_linspace.slice(1,-1);
  createAreasForImageMap(mapEl, output_linspace,"out");
  addElementToContainer(mapEl);
  return mapEl;
}

function createAreasForImageMap(mapEl, spacing, area_type="in"){
// creates areas at the specified locations; if area_type == "in", areas are created on the left
// Otherwise, they are on the right.
// mapEl: HTMLElement - HTML map element that should have the areas
// spacing: Array - array containing the y-center of the area to be added
// area_type: str - {"in", "out"}; if "in", place areas on the left of the element; else, on the right
    let nodeNum = extractIdNums(mapEl.id);
  // put area on the left
    let x_start, x_end, iarea_prefix;
    if(area_type === "in"){
        x_start = 0;
        x_end = XPUT_WIDTH;
        iarea_prefix = IAREA_IN;
    }
    else{  // put area on the right
        x_start = (NODE_WIDTH-XPUT_WIDTH).toString();
        x_end = (NODE_WIDTH).toString();
        iarea_prefix = IAREA_OUT;
    }
    // create area; add to the specified map element.
    let image_area, y_start, y_end;
    for(let idx=0; idx<spacing.length; idx++){
        image_area = document.createElement("area");
        image_area.id = formatId(iarea_prefix, nodeNum, idx);
        y_start = (spacing[idx] - XPUT_WIDTH/2).toString();
        y_end = (spacing[idx] + XPUT_WIDTH/2).toString();
        image_area.coords = x_start + "," + y_start + "," + x_end + "," + y_end;
        image_area.onmouseenter = dotHandler;
        image_area.onmouseout = removeDot;
        if(area_type !== "in") {
            image_area.onmousedown = connectorMaker;
            image_area.style.cursor = "crosshair";
        }
        mapEl.appendChild(image_area);
    }
}

function dotHandler(event){
    if(movingNode){
        return;
    }
    let area = event.target;
    let area_coords = area.coords;
    let areaId = extractIdNums(area.id);
    let areaInOut = areaId[0][0];
    let areaRect = getAreaBoundingClientRect(area);
    let mapX = areaRect.left;  // these coords are for the image map
    let mapY = areaRect.top;
    mapX += areaRect.width/2;
    mapY += areaRect.height/2;

    let xOffset;
    if(areaInOut === "in"){
    xOffset = -5;
    color = "#ee0000";
    }
    else{
    xOffset = +5;
    color = "#5500aa";
    }
    let [containerX, containerY] = getContainerOffset();
    mapX -= containerX;
    mapY -= containerY;
    let dotEl = hoverDot(mapX+xOffset, mapY, color);
    // move center of span to center of area
    let dotElRect = dotEl.getBoundingClientRect();

    dotEl.move(-dotElRect.width/2, -dotElRect.height/2);
}

function getAreaBoundingClientRect(area) {
    const img = document.querySelector(`img[usemap="#${area.parentNode.id}"]`);
    const imgRect = img.getBoundingClientRect();
    const coords = area.coords.split(',').map(Number);
    let areaRect;
    areaRect = {
        x: imgRect.x + coords[0],
        y: imgRect.y + coords[1],
        width: coords[2] - coords[0],
        height: coords[3] - coords[1],
    };
    return {
        x: areaRect.x,
        y: areaRect.y,
        width: areaRect.width,
        height: areaRect.height,
        top: areaRect.y,
        right: areaRect.x + areaRect.width,
        bottom: areaRect.y + areaRect.height,
        left: areaRect.x,
    };
}

function parseCoords(s){
// parses a coordinates string of the form x_min,y_min,x_max,y_max and returns each of these values
    let coords = [];
    let coordSplit = s.split(',');
    for(let idx=0; idx<4; idx++){
        coords[idx] = parseFloat(coordSplit[idx]);
    }
    return coords;
}

function getPscsNodeFromNum(num){
    const nodeId = formatId(NODE, num);
    return document.getElementById(nodeId);
}

function removeDot(event){
    let dotEl = document.getElementById("hoverdot");
    if(dotEl != null) {
        dotEl.style.transitionDuration = "0.5s";
        dotEl.style.transitionProperty = "opacity";
        dotEl.style.opacity = "0%";
        dotEl.addEventListener("transitionend", function removeWhenGone() {
          this.remove();
        });
    }
}

function createDot(){
// creates the span element
  let dotEl = document.createElement("span");
  dotEl.id = "hoverdot";
  dotEl.class = "hoverdot2";
  dotEl.style.height = "15px";
  dotEl.style.width = "40px";
  dotEl.style.borderRadius = "50%";
  dotEl.style.position = "absolute";
  addElementToContainer(dotEl);
  // document.body.appendChild(dotEl);
  return dotEl;
}

function hoverDot(posX, posY, color="#ff0000"){
// creates a dot at the specified position
// posX: float - X position where to place the dot
// posY: float - Y position where to place the dot
// color: string - color and transparency of dot
  let dotEl = document.getElementById("hoverdot");
  if(dotEl == null){dotEl = createDot()}
  dotEl.moveTo = moveElementToPos;
  dotEl.move = moveElement;
  dotEl.moveTo(posX,posY);
  dotEl.style.backgroundColor = color;
  dotEl.transitionDuration = "0s";
  dotEl.style.zIndex = "1";
  dotEl.style.opacity = "0.5";
  return dotEl;
}

function formatConnectorId(srcAreaId, dstAreaId=null){
  const srcIds = extractIdNums(srcAreaId);  // [0]-in, [1]-nodeNum, [2]-areaNum
  let dstIds;
  if(dstAreaId == null){
    dstIds = ["out", ".", "."];
  }
  else{
    dstIds = extractIdNums(dstAreaId);  // [0]-out, [...]
  }
  let connectorId = CONNECTOR + SEP + srcIds[1] + FIELDSEP + srcIds[2];
  connectorId += SEP + dstIds[1] + FIELDSEP + dstIds[2];
  return connectorId;
}

function connectorMaker(event){
    event.preventDefault();
    const targetArea = event.target;
    let connector = createConnector(targetArea.id);
    connector.grab(event);

}

function createConnector(srcAreaId){
  let connector = document.createElement("span");
  connector.class = CONNECTOR;
  // connector.id = formatConnectorId(srcAreaId);
  applyClassCSS(connector);

  // set connector params
  let [centerX, centerY] = getCenterOfArea(srcAreaId);

  connector.moveTo = moveElementToPos;
  connector.moveTo(centerX, centerY);
  connector.srcPos = [centerX, centerY];
  connector.dstPos = [null, null];
  connector.containerOffset = getContainerOffset();
  connector.update_pos = update_pos;
  connector.grab = grab;
  connector.del = del;
  connector.connect = connect;
  connector.disconnect = disconnect;
  connector.moveSrc = moveSrc;
  connector.moveDst = moveDst;
  connector.coMove = [];
  connector.onmousedown = function(){
    if(lastSelectedElement !== null) {
        lastSelectedElement.coMove = lastSelectedElement.coMove.filter(el => el.id !== "selectionMarker");
        document.getElementById("selectionMarker").remove();
      }
      lastSelectedElement = connector;
    let selectionMarker = makeSelectionMarker(10,10);
    alignElementCenterToOther(selectionMarker, connector);
      connector.coMove.push(selectionMarker);
  };
  addElementToContainer(connector);
  return connector;
  function grab(event){
      // move dst to mouse
      event.preventDefault();
      document.onmousemove = drag;
      document.onmouseup = release;
  }
  function drag(event){
      event.preventDefault();
      let mouseX = event.clientX - connector.containerOffset[0];
      let mouseY = event.clientY - connector.containerOffset[1];
      connector.dstPos = [mouseX, mouseY];
      connector.update_pos();
  }
  function moveSrc(moveX, moveY){
      connector.srcPos = [connector.srcPos[0] + moveX, connector.srcPos[1] + moveY];
      connector.update_pos();
  }
  function moveDst(moveX, moveY){
      connector.dstPos = [connector.dstPos[0] + moveX, connector.dstPos[1] + moveY];
      connector.update_pos();
  }
  function release(event){
      event.preventDefault();
      // on release, line should attach to nearest input area
      // get all candidate pscsNodes
      const nearestArea = getNearestInputArea(event.clientX, event.clientY);
      if(nearestArea == null){
          connector.del();
      }
      else{
          let tentativeId = formatConnectorId(srcAreaId, nearestArea.id);
          if(connectsToSelf(tentativeId)){
              connector.del();
          }
          else {
              let otherConnector = document.getElementById(tentativeId);
              if (otherConnector == null) {
                  // attach dst to the area
                  connector.id = tentativeId;
                  const [centerX, centerY] = getCenterOfArea(nearestArea);
                  // add connector to dstNode
                  connector.connect(srcAreaId, false);
                  connector.connect(nearestArea.id);
                  connector.dstPos = [centerX, centerY];
                  connector.update_pos();
              } else {
                  connector.del()
              }
          }
      }
      document.onmousemove = null;
      document.onmouseup = null;
  }
  function del(){
      // need to look at dstNode, srcNode and remove
      connector.disconnect();
      connector.remove();
      let delMarker = document.getElementById("selectionMarker");
      if(delMarker !== null){
          delMarker.remove();
      }
  }

  function disconnect(){
      // removes this connector from its attached nodes
      if(connector.id == null || connector.id.length === 0){return}
      const idNums = extractIdNums(connector.id);
      const srcNodeNum = idNums[0][0];
      const dstNodeNum = idNums[1][0];
      let srcNode = getPscsNodeFromNum(srcNodeNum);
      // remove connector ID from srcNode
      let srcIdx = srcNode.srcConnectors.indexOf(connector.id);
      srcNode.srcConnectors.splice(srcIdx, 1);
      let dstNode = getPscsNodeFromNum(dstNodeNum);
      if(dstNode != null) {
          let dstIdx = dstNode.dstConnectors.indexOf(connector.id);
          dstNode.dstConnectors.splice(dstIdx, 1);
      }
  }
  function connect(areaId, isDst=true){
      // adds the connector to the node's data; if isDst==true, adds to dstConnectors; else, to srcConnectors
      let pscsNodeNum = extractIdNums(areaId)[1][0];
      let pscsNode = getPscsNodeFromNum(pscsNodeNum);
      if(isDst){
          if(!pscsNode.dstConnectors.includes(connector.id)) {
              pscsNode.dstConnectors.push(connector.id);
          }
      }
      else{
          if(!pscsNode.srcConnectors.includes(connector.id)) {
              pscsNode.srcConnectors.push(connector.id);
          }
      }
  }

  function update_pos(){
      // move element position to match its src and dst positions
      connector.moveTo(connector.srcPos[0], connector.srcPos[1]);
      // calculate height
      let height = Math.pow(connector.dstPos[1]-connector.srcPos[1], 2) + Math.pow(connector.dstPos[0]-connector.srcPos[0],2);
      height = Math.sqrt(height);

      // calculate angle
      const radAngle = getRadAngleFromVertical(connector.dstPos[1] - connector.srcPos[1], connector.dstPos[0] - connector.srcPos[0]);
      const degAngle = radToDeg(radAngle);

      connector.style.height = height.toString() + "px";
      connector.style.transform = "rotate(" + degAngle.toString() + "deg)";
  }
}

function connectsToSelf(connectorId){
    let connectorSplit = connectorId.split(SEP);
    let firstNode = connectorSplit[1].split(FIELDSEP)[0];
    let secondNode = connectorSplit[2].split(FIELDSEP)[0];
    return firstNode === secondNode;


}


function getNearestInputArea(posX, posY){
    // Returns the image area closest to the input x,y position
    const nearestCandidates = document.elementsFromPoint(posX, posY);
    // check for imagemaps
    let keptCandidates = [];
    let mapId;
    nearestCandidates.forEach(candidate => {
        if(candidate.id.includes(NODE)){
            mapId = candidate.useMap.slice(1);
            let map = document.getElementById(mapId);
            keptCandidates.push(map);
            // keptCandidates.push(candidate.useMap.slice(1));
        }
    })
    // check which input area is closest
    let minDist = Infinity;
    let closestArea = null;
    for(let idx=0; idx<keptCandidates.length; idx++){
        // get input areas
        let areas = keptCandidates[idx].children;
        for(let areaIdx=0; areaIdx<areas.length; areaIdx++){
            // ensure that this is an input area
            let distance;
            if(areas[areaIdx].id.includes("in")){
                distance = computeDistanceToArea(areas[areaIdx], posX, posY);
                if(distance < minDist){
                    minDist = distance;
                    closestArea = areas[areaIdx];
                }
            }
        }
    }
    return closestArea;
}

function computeDistanceToArea(area, posX, posY){
    let [centerX, centerY] = getCenterOfArea(area, false);
    let distance = Math.pow(posY - centerY,2) + Math.pow(posX-centerX,2);
    return Math.sqrt(distance);
}

function getContainer(containerId = "nodeContainer"){
    return document.getElementById(containerId);
}

function getContainerOffset(containerId = "nodeContainer"){
    let container = document.getElementById(containerId);
    let containerRect = container.getBoundingClientRect();
    return [containerRect.left, containerRect.top];
}

function getRadAngleFromVertical(deltaY, deltaX){
    // Gets the angle of rotation for a displacement relative to the vertical axis.
    let atan = Math.atan(deltaY / deltaX)
    if (deltaX < 0){
        atan += Math.PI/2;
    }
    else {
        atan -= Math.PI/2;
    }
    return atan;
}

function radToDeg(rad) {
    return rad * 180 / Math.PI;
}

function getCenterOfArea(area, containerOffset = true){
    // get coords
    if(typeof(area) == 'string'){
        area = document.getElementById(area);
    }
    let areaRect = getAreaBoundingClientRect(area);
    let mapX = areaRect.left + areaRect.width/2;
    let mapY = areaRect.top + areaRect.height/2;
    if(containerOffset){
        const [contX, contY] = getContainerOffset()
        return [mapX-contX, mapY-contY]
    }
    else{return [mapX, mapY]}
}



function addElementToContainer(element, containerId="nodeContainer"){
// central function for adding new elements to the page
  const container = document.getElementById(containerId);
  container.appendChild(element);
}

function parseParams(paramString){
  // params: string - JSON string of form {paramName: [paramType, paramValue]}
    if(typeof(paramString) === "string") {
        paramString = paramString.replaceAll("'", '"');
        paramString = paramString.replaceAll("None", null);
        return JSON.parse(paramString);
    }
    else{
        return paramString;
    }

}

function linspace(startValue, endValue, numValues=2){
  if(numValues === 1){
    return (endValue-startValue)/2;
  }
  let spacing = (endValue - startValue) / (numValues-1);
  let lin = [];
  for(let idx=0; idx<numValues; idx++){
    lin.push(spacing*idx);
  }
  return lin;
}

function applyClassCSS(element, cls = "") {
// Fetches the CSS for the class of the element and applies it. If class is specified, apply the
//CSS for that class instead.
//element : HTMLElement
//    Elements to which to apply the CSS.
//class : String
//    Optional. Class for which to fetch the CSS rule; takes precedence over the element's class.
//
// Determine which class to check
let class_to_check = element.class;
if (cls.length !== 0) {
  class_to_check = cls;
}
// Find the matching rule
for (const sheet of document.styleSheets) {
    for (const r of sheet.cssRules){
        let splitText = r.selectorText.replace(" ","").split(",");
        if (splitText.includes("." + class_to_check)){
            element.style = r.style.cssText;
        }
    }
  }
}

async function exportPipeline(inputNameElement, pipelineDestination, panel) {
// This function pulls together all nodes, their properties, and connections
    let pipelineName = inputNameElement.value.toString()
    if(pipelineName.length === 0){
        inputNameElement.setCustomValidity("The pipeline must have a name.");
        inputNameElement.reportValidity()
        return;
    }
    let isDestProject = pipelineDestination.startsWith('project:');
    let dest = pipelineDestination.substring(pipelineDestination.indexOf(':')+1);
    let node_export = [];
    for(let i=0; i<nodeIds.length; i++){
        let el = document.getElementById(nodeIds[i]);
        node_export.push(el);
    }
    let pipelineSummary = {};
    pipelineSummary['name'] = pipelineName;
    pipelineSummary['nodes'] = node_export;
    pipelineSummary['isDestProject'] = isDestProject;
    pipelineSummary['saveId'] = dest;
    await fetch("/pipeline", {
                        method: "POST",
                        headers: {
                            "Accept": "application/json",
                            "Content-Type": "application/json"
                        },
                        body: JSON.stringify(pipelineSummary)
                    });
    panel.remove();
    alert("Pipeline saved! You can execute the pipeline from the pipeline page.");
}

function openSavePanel(projDestsJSON) {
    projDestsJSON = projDestsJSON.replaceAll("'", '"');
    projDestsJSON = projDestsJSON.replaceAll("None", null);
    let projDests = JSON.parse(projDestsJSON);

    const panel = document.createElement("div");
    panel.class = EXPORTPANEL;
    panel.id = formatId(EXPORTPANEL, "save");
    // add below container
    let containerRect = getContainerOffset();
    panel.style.left = (containerRect[0]).toString() + "px";
    applyClassCSS(panel);

    let tbl = document.createElement('table');
    tbl.style.width = '100px';
    tbl.style.border = '5px solid black';

    // Need: pipeline name, associated project
    let tr = tbl.insertRow();
    let td = tr.insertCell();
    let inpName = document.createElement("INPUT");
    inpName.placeholder = "Name of the pipeline";
    td.appendChild(inpName);

    // Select menu for project destinations
    td = tr.insertCell();
    let drop = document.createElement("select");
    drop.style.maxWidth = '100px';
    let option = document.createElement("option");
    option.innerHTML = "Destination";
    option.disabled = true;
    option.selected = true;
    drop.appendChild(option);

    // go through all projects with which the user is associated
    for (let proj in projDests){
        option = document.createElement("option");
        let projDisplay = proj;
        if (proj.length > 20){
            projDisplay = proj.substring(0,17) + "...";
        }
        option.innerHTML = "Project: " + projDisplay;
        option.value = 'project:' + projDests[proj];
        drop.appendChild(option);
    }
    td.appendChild(drop);

    // Save button
    td = tr.insertCell();
    let btn = document.createElement("button");
    btn.innerHTML = "Save";
    btn.onclick = function() {exportPipeline(inpName, drop.value.toString(), panel);}
    td.appendChild(btn);

    // Close button
    td = tr.insertCell();
    btn = document.createElement("button");
    btn.innerHTML = "&times;";
    btn.onclick = function() {panel.remove();}
    td.appendChild(btn)

    panel.appendChild(tbl);
    document.body.appendChild(panel);
}

async function loadAnalysis(){
    let analysis_selection = document.getElementById("analyses");
    let id_analysis = analysis_selection.value;
    return loadAnalysisFromId(id_analysis);
}


async function loadAnalysisFromId(id_analysis){
// Loads analysis data from the DB based on the input ID.
    let analysis_spec = {};
    analysis_spec["loadAnalysis"] = id_analysis;
    let response = await fetch(window.location.href + "/load_analysis", {
        method: "POST",
        headers: {
            "Accept": "application/json",
            "Content-Type": "application/json"},
        body: JSON.stringify(analysis_spec)});
    let data = await response.json();
    let nodes = data["nodes"];
    clearNodes();
    let nodeObject = {};
    let nodesWithoutDepth = [];
    for(let idx=0; idx<nodes.length; idx++){
        let node = nodes[idx];
        let load_params = {};
        let paramName, paramType, paramValue;
        for(let param_idx=0; param_idx<node["params"].length; param_idx++){
            paramName = node["params"][param_idx];
            paramType = node["paramsTypes"][paramName];
            paramValue = node["paramsValues"][paramName];
            load_params[paramName] = [paramType, paramValue];
        }
        let pkgName = node["module"].split(".")[0]
        let pkg = selectPackage(pkgName, nodeInfo);
        let psNode = getNodeDataFromModule(node["module"], node["procName"], pkg["modules"]);

        let nodeNum = extractIdNums(node["nodeId"]);

        let loadedNode = createPscsNode(nodeNum, node["img"], psNode);
        for(let paramName in load_params){
            loadedNode.paramsValues[paramName] = load_params[paramName][1];
        }
        if(node.pscsType === "input"){
            loadedNode.depth = 0;
        }
        else{
            loadedNode.depth = -1;
            nodesWithoutDepth.push(loadedNode);
        }
        loadedNode.id = node.nodeId;
        loadedNode.dstConnectors = node.dstConnectors;
        loadedNode.srcConnectors = node.srcConnectors;
        nodeNum = extractIdNums(node.nodeId)[0];
        nodeObject[nodeNum] = loadedNode;
    }

    // make connectors
    // get container for reference
    const container = document.getElementById("nodeContainer");
    const containerRect = container.getBoundingClientRect();
    let currentNode, currentNodeNum, srcNodeLeft, srcNodeTop, srcString;
    let connectorId, parsedConnectorId, srcNodeNum, srcAreaNum, dstNodeNum, dstAreaNum;
    let srcAreaId, srcCenter, dstAreaId, dstCenter, connector;
    for(let nodeIdx in nodeObject){
        // get the node being considered
        currentNode = nodeObject[nodeIdx];
        currentNodeNum = extractIdNums(currentNode.id)[0];
        const currentNodeRect = currentNode.getBoundingClientRect();
        // go through src
        srcNodeLeft = currentNodeRect.left - containerRect.left;
        srcNodeTop = currentNodeRect.top - containerRect.top;
        for(let srcIdx in currentNode.srcConnectors){
            srcString = currentNode.srcConnectors[srcIdx];
            connectorId = currentNode.srcConnectors[srcIdx];
            parsedConnectorId = extractIdNums(connectorId);
            srcNodeNum = parsedConnectorId[0][0];
            srcAreaNum = parsedConnectorId[0][1];
            dstNodeNum = parsedConnectorId[1][0];
            dstAreaNum = parsedConnectorId[1][1];
            srcAreaId = formatId(IAREA_OUT, srcNodeNum, srcAreaNum);
            srcCenter = getCenterOfArea(document.getElementById(srcAreaId), true);
            dstAreaId = formatId(IAREA_IN, dstNodeNum, dstAreaNum);
            dstCenter = getCenterOfArea(document.getElementById(dstAreaId), true);
            connector = createConnector(srcAreaId);
            connector.id = formatConnectorId(srcAreaId, dstAreaId);
            connector.connect(srcAreaId, false);
            connector.connect(dstAreaId);
            connector.srcPos = srcCenter;
            connector.dstPos = [dstCenter[0], dstCenter[1]];
            connector.update_pos();
            connector.srcNode = formatId(NODE, srcNodeNum);
            connector.dstNode = formatId(NODE, dstNodeNum);
            document.onmousemove = null;
            document.onmouseup = null;
        }
    }
    // Below this is aesthetics.
    // determine depth of each node
    determineDepths();
    let nodeList = document.querySelectorAll("[id^=" + NODE + "]");
    let maxDepth = 0;
    nodeList.forEach(el => {if(el.depth > maxDepth){maxDepth=el.depth}});
    // determine x position of each
    let xposByDepth = {};
    for(let idx=0; idx<maxDepth+2; idx++){
        xposByDepth[idx] = 0;
    }
    let node;
    for(let n in nodeObject){
        node = nodeObject[n];
        node.xpos = xposByDepth[node.depth];
        xposByDepth[node.depth]+=1;
    }

    // move nodes to position
    let xpos, ypos;
    for(let idx in nodeObject){
        node = nodeObject[idx];
        ypos = node.xpos * HORIZONTAL_SPACING + Math.random() * HORIZONTAL_JIGGLE;
        xpos = node.depth * VERTICAL_SPACING + Math.random() * VERTICAL_JIGGLE;
        node.moveAll(xpos, ypos);
    }

}

function determineDepths(){
    let nodeList = document.querySelectorAll("[id^=" + NODE + "]");
    let inputNodes = [];
    nodeList.forEach(el => {
        if(el.pscsType === "input"){
            inputNodes.push(el);
        }
    });
    // Check the nodes following the input nodes
    let nodesToCheck = [];
    inputNodes.forEach(el => {
            nodesToCheck = nodesToCheck.concat(getNextNodes(el));
        });

    do{
        let nextCheck = [];
        for(let checkIdx=0; checkIdx<nodesToCheck.length; checkIdx++){
            let nodeToCheck = nodesToCheck[checkIdx];
            let prevDepth = maxDepthOfPrevNodes(nodeToCheck);
            if(prevDepth === -1){
                // instead, find the previous node that doesn't have depth, then check it
                let prevNodes = getPrevNodes(nodeToCheck);
                nextCheck = nextCheck.concat(getDepthless(prevNodes));
            }
            else{
                nodeToCheck.depth = prevDepth + 1;
                // check next ones
                let nextNodes = getNextNodes(nodeToCheck);
                nextCheck = nextCheck.concat(getDepthless(nextNodes));
            }
        }
        nodesToCheck = nextCheck;
    } while(nodesToCheck.length > 0);
}

function hasDepth(nodeToCheck){
    return !(nodeToCheck.depth === -1 || nodeToCheck.depth === undefined)
}

function getDepthless(nodeList){
    // checks each element of nodeList to verify that they have .depth >= 0
    let depthless = [];
    nodeList.forEach(el =>{
        if(!hasDepth(el)){
            depthless.push(el);
        }
    });
    return depthless;
}

function maxDepthOfPrevNodes(pscsNode){
    // Checks the depth of previous nodes; returns -1 if .depth is undefined or -1
    // iterate through dstConnector
    let maxDepth = -1;
    let prevNodes = getPrevNodes(pscsNode);
    for(let idx=0; idx<prevNodes.length; idx++){
        let prevNode = prevNodes[idx];
        if(!hasDepth(prevNode)){
            return -1
        }
        else{
            maxDepth = Math.max(maxDepth, prevNode.depth);
        }
    }
    return maxDepth
}

function getNextNodes(pscsNode){
    // Returns an array of nodes that come after the input node.
    let nextNodes = []
    pscsNode.srcConnectors.forEach(conn => {
        let nextNum = extractIdNums(conn)[1][0];  // note; if you want to track port, this is what you need to change
        nextNodes.push(getPscsNodeFromNum(nextNum));
    })
    return nextNodes;
}

function getPrevNodes(pscsNode){
    let prevNodes = []
    pscsNode.dstConnectors.forEach(conn => {
        let prevNum = extractIdNums(conn)[0][0];
        prevNodes.push(getPscsNodeFromNum(prevNum));
    })
    return prevNodes;
}

function clearNodes() {
    // removes all nodes from the page
    let nodeId, node;
    for (let idx = nodeIds.length - 1; idx >= 0; idx--) {
        nodeId = nodeIds[idx];
        node = document.getElementById(nodeId);
        if(node != null) {
            node.del();
        }
    }
    nodeIds = [];
}

async function getNodeInfo(){
    // fetches the node information from the server
    let response = await fetch("/pipeline/fetch_nodes", {
        method: "POST",
        headers: {
            "Accept": "application/json",
            "Content-Type": "application/json"
        }
    })
    nodeInfo = await response.json()
    nodeInfo = nodeInfo["packages"];
    document.dispatchEvent(new Event('nodesLoaded'));
    return nodeInfo;
}

function makeSidebarNode(parentDiv, nodeInfo, depth){
    // make HTML element
    let p = document.createElement("p");
    p.title = JSON.stringify(nodeInfo["name"]);
    p.style.paddingLeft = String(depth*PADDING_LEVEL) + "px";
    let nodeType = nodeInfo["type"];
    const imgPath = "/static/nodes/" + nodeType + ".png";
    p.addEventListener('click',
                        function () {createPscsNode(null, imgPath, nodeInfo);});
    p.textContent = nodeInfo["name"];
    p.style.fontSize = getFontSize(depth);
    p.style.cursor = "pointer";
    let img = document.createElement("img");
    img.src = imgPath;
    img.width = 20;
    img.height = 20;
    p.appendChild(img);
    parentDiv.appendChild(p);
}
let unwrappedNodes = [];
function makeNestedSidebarModule(parentDiv, moduleInfo, depth){
    // adds the module to the sidebar; if the module has modules, adds those, too
    let top_div = addCollapsibleSection(parentDiv, moduleInfo["name"], depth);
    for(let n of moduleInfo["nodes"]){
        makeSidebarNode(top_div, n, depth+1);
        unwrappedNodes.push(n);
    }
    for(let m of moduleInfo["modules"]) {
        makeNestedSidebarModule(top_div, m, depth+1);
    }
}

async function populateSidebar(sidebarId) {
    // Adds node info to the sidebar
    for (let pkg of nodeInfo) {
        let pkg_name = pkg["display_name"];
        let modules = pkg["modules"];
        let top_div = document.getElementById(sidebarId);
        let depth = 0
        top_div = addCollapsibleSection(top_div, pkg_name, depth);
        modules = modules["modules"];
        for (let m of modules) {
            makeNestedSidebarModule(top_div, m, depth+1);
        }
    }
}

function getFontSize(depth){
    return String(Math.max(START_FONTSIZE-DECREASE_STEP_FONTSIZE*depth, MIN_FONTSIZE)) + "px";
}

function addCollapsibleSection(parentDiv, title="", depth){
    // Adds a collapsible subsection to the div specified by divId
    let section_title = document.createElement("p");
    section_title.textContent = title;
    section_title.style.fontSize = getFontSize(depth);
    section_title.addEventListener("click", toggleDisplay);
    section_title.style.paddingLeft = String(depth*PADDING_LEVEL) + "px";
    section_title.style.cursor = "pointer";
    parentDiv.appendChild(section_title);
    let collap = document.createElement("div");
    collap.classList.add("sidebar-section");
    collap.style.paddingLeft = String(depth*PADDING_LEVEL) + "px";
    parentDiv.appendChild(collap);
    section_title.click();  // collapse section
    return collap;
}

function toggleDisplay(e){
    // Fires on event; toggles display of sibling
    let el = e.target;
    let sib = el.nextElementSibling;
    if(sib.style.display === "none"){
        sib.style.display = "block";
        if(el.textContent[0] === "+"){
            el.textContent = "-" + el.textContent.substring(1,);
        }
        else{
            el.textContent = "-" + el.textContent;
        }
    }
    else{
        sib.style.display = "none";
        if(el.textContent[0] === "-"){
            el.textContent = "+" + el.textContent.substring(1,);
        }
        else{
            el.textContent = "+" + el.textContent;
        }
    }
}

function getNodeType(pscsNode){
    // Determines whether the node is input, output, simo, or mimo
    if(pscsNode["num_inputs"] === 0){return "input"}
    else if(pscsNode["num_outputs"] === 0){return "output"}
    else if(pscsNode["num_inputs"] === 1){return "simo"}
    else{return "mimo"}
}

function setTextContent(elementID, text){
    document.getElementById(elementID).textContent = text;
}

function startImport(projectListID, importButtonID){
  let projectList = document.getElementById(projectListID);
  let importButton = document.getElementById(importButtonID);
  projectList.style.display = "block";
  importButton.style.display = "block";
}

function validatePipeline(labelId){
    // Verifies that all nodes have the appropriate parameters set.
    // First remove any existing markers
    let markers = document.getElementsByClassName("invalidMarker");
    while(markers.length > 0){
        markers[0].remove();
    }
    let allNodes = document.getElementsByClassName(NODE);
    let pipelineValid = true;
    // Resolve parameters for all nodes.
    for(let node of allNodes){node.resolveInteractions()}
    for(let node of allNodes){node.cumulativeEffect()}

    for(let i=0; i<allNodes.length; i++){
        let [nodeIsValid, invalidReasons] = validateNode(allNodes[i]);
        if(!nodeIsValid){
            pipelineValid = false;
            allNodes[i].setInvalid(invalidReasons.join("\n"));
        }
    }
    if(!pipelineValid){
        document.getElementById(labelId).style.display = "block";
    }
    else{
        document.getElementById(labelId).style.display = "none";
    }

}

function validateNode(nodeEl){
    // Verifies that a node has the appropriate parameters set, isn't hanging, etc.
    let undefinedParams = getUndefinedRequiredParameters(nodeEl);
    let inputsConnected = areInputsConnected(nodeEl);
    let outputsConnected = areOutputsConnected(nodeEl);
    let uniquePortConnections = doPortsHaveOneConnection(nodeEl);
    let requirementsMet = nodeEl.evaluateReady();
    let isValid = (undefinedParams.length === 0) && inputsConnected && outputsConnected && uniquePortConnections && requirementsMet;
    let invalidReasons = [];
    if(!isValid) {
        if (undefinedParams.length !== 0) {
            invalidReasons.push("Undefined parameters: " + undefinedParams.join(", ") + ".");
        }
        if (!inputsConnected) {
            invalidReasons.push("Not all inputs are connected.");
        }
        if (!outputsConnected) {
            invalidReasons.push("Not all outputs are connected.");
        }
        if (!uniquePortConnections) {
            invalidReasons.push("Input ports have multiple connections.");
        }
        if (!requirementsMet) {
            invalidReasons.push(prepareMissingReqMessage(nodeEl.unmetRequirements[0]));
        }
    }
    return [isValid, invalidReasons]
}

function prepareMissingReqMessage(missingReqs){
    let msg = "The following fields are required by the node but are not produced by the pipeline. Either nodes " +
        "are missing before this one or the fields must be present in the input data.\n"
    for(let key in missingReqs){
        if(missingReqs[key].length > 0){
            msg += `${key}: ${missingReqs[key]}\n`
        }
    }
    let recNodes = getNRecommendedEffectNodes(missingReqs, 3);
    if(recNodes.length > 0) {
        msg += "The following nodes might be helpful:\n"
        for (let rn of recNodes) {
            msg += `Under ${rn["module"]}: ${rn["name"]}\n`
        }
    }
    return msg;
}

function getNRecommendedEffectNodes(reqs, n=3){
    let effNodes = findEffectNode(reqs);
    let recNodes = [];
    let count = 0;
    let idx = 0;
    let exnodes = effNodes[idx];
    while(count < n && idx < effNodes.length){
        recNodes = recNodes.concat(exnodes.slice(0, n-count));
        count = recNodes.length;
        idx += 1;
        exnodes = effNodes[idx];
    }
    return recNodes;
}


function findEffectNode(reqs){
    let numMissingEffects = getInteractionCount(reqs);
    let recommendedNodes = [];
    for(let i=0; i<numMissingEffects; i++){recommendedNodes.push([]);}
    for(let node of unwrappedNodes){
        if(node.effects === undefined || node.effects === null || node.effects.length === 0){
            continue
        }
        let overlap = getInteractionOverlapCount(reqs, node.effects[0]);
        if(overlap === 0){
            continue
        }
        recommendedNodes[numMissingEffects-overlap].push(node);
    }
    return recommendedNodes;
}

function getInteractionCount(interactions){
    let count = 0;
    for(let k in interactions){
        count += interactions[k].length;
    }
    return count;
}

function getInteractionOverlapCount(inter0, inter1){
    let overlapCounter = 0;
    for(let k in inter0){
        for(let val0 of inter0[k]) {
            if (inter1[k].includes(val0)) {
                overlapCounter += 1;
            }
        }
    }
    return overlapCounter;
}


function getUndefinedRequiredParameters(nodeEl){
    let undefinedParams = [];
    let requiredParams = nodeEl.required_parameters;
    for(let i = 0; i < requiredParams.length; i++){
        if(nodeEl.paramsValues[requiredParams[i]] === null){
            undefinedParams.push(requiredParams[i]);
        }
    }
    return undefinedParams
}

function areInputsConnected(nodeEl){
    // For each port, there should be an entry in .dstConnectors
    let portsWithConnections = new Set();
    let idInfo = "";
    let this_port = "";
    for(let i=0; i<nodeEl.dstConnectors.length; i++){
        this_port = extractIdNums(nodeEl.dstConnectors[i])[1][1];
        // not checking if port has multiple connections in this function
        portsWithConnections.add(this_port);
    }
    return portsWithConnections.size === nodeEl.num_inputs;
}

function doPortsHaveOneConnection(nodeEl){
    let inputPortConnections = new Set();
    let this_port = "";
    for(let i= 0; i<nodeEl.dstConnectors.length; i++){
        this_port = extractIdNums(nodeEl.dstConnectors[i])[1][1];
        if(inputPortConnections.has(this_port)){
            return false;
        }
        else {
            inputPortConnections.add(this_port);
        }
    }
    return true;
}

function areOutputsConnected(nodeEl){
    // For each port, there should be an entry in .dstConnectors
    let portsWithConnections = new Set();
    let idInfo = "";
    let this_port = "";
    for(let i=0; i<nodeEl.srcConnectors.length; i++){
        this_port = extractIdNums(nodeEl.id)[0][1];
        // not checking if port has multiple connections in this function
        portsWithConnections.add(this_port);
    }
    return portsWithConnections.size === nodeEl.num_outputs;
}

async function importPipeline(analysisElemId, projectListElemID){
    let analysisEl = document.getElementById(analysisElemId);
    let projectEl = document.getElementById(projectListElemID);
    let project = projectEl.options[projectEl.selectedIndex].value;
    let analysis = analysisEl.options[analysisEl.selectedIndex].value;
    let projectName = projectEl.options[projectEl.selectedIndex].text;
    let analysisName = analysisEl.options[analysisEl.selectedIndex].text;
    let importData = {"id_project": project, "id_analysis": analysis};
    let response = await fetch("./public_pipeline_import", {
                                            method: "POST",
                                            headers: {
                                                "Accept": "application/json",
                                                "Content-Type": "application/json"
                                            },
                                            body: JSON.stringify(importData)
                                            });
    let resp = await response.json();
    if(resp["success"] === 1){
        alert("Pipeline \"" + analysisName + "\" has been imported to project \"" + projectName + "\"");
    }
    else{
        alert("Pipeline not imported: " + response["message"]);
    }
}

async function requestCheckedFiles(checkClass){
    let checks = document.getElementsByClassName(checkClass);
    let requestedFiles = [];
    for(let i = 0; i < checks.length; i++) {
        if(checks[i].checked) {
            requestedFiles.push(checks[i].id);
        }
    }
    await requestFiles(requestedFiles);
}

async function requestFiles(requestedFiles) {
    // For elements of the class `checkClass` that have been checked, gets their ID and sends a request to the server
    // to download it.
    await fetch(window.location.pathname + "/file_request", {
                                            method: "POST",
                                            headers: {
                                                "Accept": "application/json",
                                                "Content-Type": "application/json"
                                            },
                                            body: JSON.stringify(requestedFiles)})
        .then(response => {
            let contentDisp = response.headers.get("Content-Disposition");
            let warningEl = document.getElementById("rateWarning");
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
        });
}

function toggleVisibility(elId){
    el = document.getElementById(elId);
    el.style.display = (el.style.display === "none") ? "block" : "none";
    return
}

function toggleIcon(elId){
    let el = document.getElementById(elId);
    if(el.innerText.startsWith("-")){
        el.innerText = el.innerText.replace("-", "+");
    }
    else if(el.innerText.startsWith("+")){
        el.innerText = el.innerText.replace("+", "-");
    }
    return
}


function makeSelectionMarker(height=null, width=null){
    let marker = document.createElement("span");
    marker.classList.add("selectionMarker");
    if(height !== null){
        marker.style.height = height + "px";
    }
    if(width !== null){
        marker.style.width = width + "px";
    }
    marker.id = "selectionMarker";
    marker.move = moveElement;
    marker.moveTo = moveElementToPos
    marker.ondblclick = lastSelectedElement.ondblclick;
    addElementToContainer(marker);
    return marker
}

function alignElementCenterToOther(el0, el1, offset=null){
    let srcRec = el0.getBoundingClientRect();
    let targetRec = el1.getBoundingClientRect();
    let targetCenterX = targetRec.left + targetRec.width/2;
    let targetCenterY = targetRec.top + targetRec.height/2;
    const [containerX, containerY] = getContainerOffset();
    if(offset === null){offset = [0,0]}
    el0.style.left = (targetCenterX - srcRec.width/2-containerX+offset[0]) + "px";
    el0.style.top = (targetCenterY - srcRec.height/2-containerY+offset[1]) + "px";
}

function delMarker(){
    let delMarker = document.getElementById("selectionMarker");
    if(delMarker !== null){
        delMarker.remove();
    }
}