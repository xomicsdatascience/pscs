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

const VERTICAL_POSITION_INITIAL = 200;  // initial vert position for new nodes (px)
const HORIZONTAL_POSITION_INITIAL = 200;  // initial horizontal position for new nodes (px)
const VERTICAL_MARGIN = 30;  // margin between new element and previous ones (px)
const HORIZONTAL_MARGIN = 30;  // margin between new element and previous ones (px)
const VERTICAL_SPACING = 100;
const HORIZONTAL_SPACING = 100;
const VERTICAL_JIGGLE = 10;
const HORIZONTAL_JIGGLE = 15;

// document listeners
document.addEventListener("keydown", kbshortcut);
function kbshortcut(event){
    if(event.keyCode === DELKEY) {
        // in case user is deleting something in a textbox:
        if(!(document.activeElement instanceof HTMLInputElement)) {
            if (lastSelectedElement != null) {
                lastSelectedElement.del();
                lastSelectedElement = null;
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

function createPscsNode(idNum, processName, module, params, pscsType, img){
  // idNum: int - optional; if defined, will attempt to make a node with this id: NODE-idNum
  // processName: string - name of the function; used for initial labelText
  // module: string - name of the Python module that the process came from
  // params: string - JSON string containing process parameters of form: paramName: [paramType, paramValue]
  // pscsType: string - type of node: {"input", "output", "simo", "mimo"}
  // img: string - path to the image to use for this node
    // Get Id for node
    let pageEl = document.createElement("img");
    pageEl.class = NODE;
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

    nodeIds.push(pageEl.id);  // add id to the list -- should probably remove this and use querySelectorAll instead
    let params_dict = parseParams(params);
    pageEl.params = [];
    pageEl.paramsValues = {};
    pageEl.paramsTypes = {};
    // store param names, values into pageEl
    for(const p in params_dict){
        pageEl.params.push(p);
        pageEl.paramsTypes[p] = params_dict[p][0];
        pageEl.paramsValues[p] = params_dict[p][1];
    }
    pageEl.title = pageEl.params.toString();
    pageEl.pscsType = pscsType;
    pageEl.img = img; // for later loading
    pageEl.src = img;
    pageEl.srcConnectors = [];
    pageEl.dstConnectors = [];

    // Create label
    pageEl.labelText = pageEl.procName;  // user can set label name, but not procName

    // Create image map
    let [numInput, numOutput] = getNumInputOutputs(pageEl);
    if(pscsType === "input"){numInput = 0}
    else if(pscsType === "output"){numOutput=0}
    pageEl.imagemap = createImageMap(pageEl.id, numInput, numOutput);
    pageEl.useMap = "#" + pageEl.imagemap.id;
    // add element to page
    pageEl.style.zIndex = "2";
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

    // define methods
    let startX = 0;
    let startY = 0;
    function grab(event){
    // signals that the element has been clicked on and should follow the cursor
      event.preventDefault();
      event.target.style.cursor = "grabbing";
      lastSelectedElement = event.target;  // for deletion
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
        let container = document.getElementById("nodeContainer");
        let containerOffset = getContainerOffset();
        let pageElRect = pageEl.getBoundingClientRect();
        let minWidth = pageElRect.right - containerOffset[0]+IAREA_MARGIN;
        let minHeight = pageElRect.bottom - containerOffset[1]+IAREA_MARGIN*2;
        container.style.minWidth = (minWidth).toString() + "px";
        container.style.minHeight = (minHeight).toString() + "px";
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
        lastSelected = null;  // prevent character deletion from user hitting 'delete'
        const panel = document.createElement("div");
        panel.class = PARAMPANEL;
        panel.id = formatId(PARAMPANEL, extractIdNums(pageEl.id)[0]);
        panel.nodeId = pageEl.id;
        applyClassCSS(panel);
        // Create table for user input
        let tbl = document.createElement('table');
        tbl.style.width = '100px';
        tbl.style.border = '5px solid black';
        tbl.style.borderCollapse = "separate";
        tbl.style.borderSpacing = "2px";
        // Name for display / input selection later
        let tr = tbl.insertRow();
        let td = tr.insertCell();
        td.appendChild(document.createTextNode('Node label'));
        td = tr.insertCell();
        let inp = document.createElement('input');
        inp.defaultValue = pageEl.labelText;
        inp.id = formatId(NODENAME, extractIdNums(pageEl.id)[0]);
        td.appendChild(inp);
        for (let idx=0; idx<pageEl.params.length; idx++) {
            // Get labels / default values, if any
            let parName = pageEl.params[idx];
            let parVal = pageEl.paramsValues[parName];
            let parType = pageEl.paramsTypes[parName];  // todo: use for validation
            // Create row; one per parameter
            tr = tbl.insertRow();
            td = tr.insertCell();
            td.appendChild(document.createTextNode(parName));
            td = tr.insertCell();
            inp = document.createElement('INPUT');
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
        // btn.onclick = closePanel;
        btn.onclick = function(){closePanel(btn.panelId)};
        btn.panelId = panel.id;
        btn.innerHTML = "&times; Cancel"
        td.appendChild(btn);
        td = tr.insertCell();
        btn = document.createElement("button");
        btn.onclick = saveParams;
        btn.nodeId = pageEl.id;
        btn.panelId = panel.id;
        btn.innerHTML = "&#10003; Save"
        td.appendChild(btn);

        panel.appendChild(tbl);
        addElementToContainer(panel);
        lastOpenPanel = panel;
    }
    function closePanel(panelId){
        // close panel; ignore any modifications
        // let btn = event.target;
        // let panel = document.getElementById(btn.panelId);
        // panel.remove();
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
    return pageEl;
}

function getNumInputOutputs(nodeEl){
    // gets the number of input/output ports
    return [1,1]  // placeholder
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
    // attach methods
    labelEl.move = moveElement;
    labelEl.del = del;
    labelEl.update = update;
    addElementToContainer(labelEl);
    return labelEl;

    function del(){
        labelEl.remove();
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
    let areaRect = area.getBoundingClientRect();
    let mapX = areaRect.left;  // these coords are for the image map
    let mapY = areaRect.top;
    area_coords = parseCoords(area_coords);
    mapX += (area_coords[0] + area_coords[2])/2;  // these move the top-left of the span to the center of the area
    mapY += (area_coords[1] + area_coords[3])/2;
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
  connector.onmousedown = function(){lastSelectedElement = connector};
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
          let otherConnector = document.getElementById(tentativeId);
          if(otherConnector == null) {
              // attach dst to the area
              connector.id = tentativeId;
              const [centerX, centerY] = getCenterOfArea(nearestArea);
              // add connector to dstNode
              connector.connect(srcAreaId, false);
              connector.connect(nearestArea.id);
              connector.dstPos = [centerX, centerY];
              connector.update_pos();
          }
          else{connector.del()}
      }
      document.onmousemove = null;
      document.onmouseup = null;
  }
  function del(){
      // need to look at dstNode, srcNode and remove
      connector.disconnect();
      connector.remove();
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
    let areaRect = area.getBoundingClientRect();
    let area_coords = parseCoords(area.coords);
    let mapX = areaRect.left;
    let mapY = areaRect.top;
    mapX += (area_coords[2] + area_coords[0])/2;
    mapY += (area_coords[3] + area_coords[1])/2;
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
    let analysis_spec = {};
    analysis_spec["loadAnalysis"] = id_analysis;
    let response = await fetch(window.location.href, {
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
        let nodeNum = extractIdNums(node["nodeId"]);
        let loadedNode = createPscsNode(nodeNum, node.labelText, node.module, load_params, node.pscsType, node.img);
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
            srcCenter = getCenterOfArea(document.getElementById(srcAreaId), false);
            dstAreaId = formatId(IAREA_IN, dstNodeNum, dstAreaNum);
            dstCenter = getCenterOfArea(document.getElementById(dstAreaId), false);
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
    let maxDepth = -1;
    let doBreak = false;
    while(nodesWithoutDepth.length > 0){
        let spliceIdxList = [];
        for(let spliceIdx = 0; spliceIdx < nodesWithoutDepth.length; spliceIdx++){
            let node = nodesWithoutDepth[spliceIdx];
            // examine previous nodes; if they all have depth, assign max depth+1 to this node
            maxDepth = -1;
            doBreak = false;
            for(let previousNodeIdx = 0; previousNodeIdx < node.dstConnectors.length; previousNodeIdx++){
                // get num of previous node
                let prevNum = extractIdNums(node.dstConnectors[previousNodeIdx])[0][0];
                let prevNode = nodeObject[prevNum];
                if(prevNode.depth === -1){
                    doBreak = true;
                    break;
                }
                // has depth; compare to maxdepth
                maxDepth = Math.max(maxDepth, prevNode.depth);
            }
            if(doBreak){
                break;
            }
            node.depth = maxDepth + 1;
            spliceIdxList.push(spliceIdx);
        }
        spliceIdxList.reverse();  // going backwards to preserve idx validity as we remove entries
        spliceIdxList.forEach(idx => {
            nodesWithoutDepth.splice(idx, 1);
        });
    }

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
}
