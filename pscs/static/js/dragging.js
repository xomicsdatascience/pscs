// This JS file contains code for dragging, stretching, and reorienting
// HTML elements.

// Global variables
var lastSelected = null;  // currently-selected HTML element
var nodeCount = 0;  // used for creating new node & map IDs
var nodeIds = new Array();
var lastOpenPanel;
// constants
const VERTICAL_POSITION_INITIAL = 200;
const HORIZONTAL_POSITION_INITIAL = 400;
const VERTICAL_SPACING = 100;
const HORIZONTAL_SPACING = 100;
const VERTICAL_JIGGLE = 10;
const HORIZONTAL_JIGGLE = 15;
// ID substrings
const SEPSTR = '-';
const DRAGSTR = "draggable";
const LINESTR = "connector";
const IMAPSTR = "imagemap";
const AREASTR = "area";
const NODECLASS = "dragged"
const TEXTCLASS = "dragtext";
const PANELCLASS = "panel";
const INPUTCLASS = "input";
const NAMECLASS = "name";
// Facilitate code reading
const DELKEY = 46;

// For removing elements.
document.addEventListener("keydown", delNode, false);
function delNode(event) {
    if (event.keyCode == DELKEY) {
      if (lastSelected == null){
        return;
      }
      lastSelected.del();
      lastSelected = null;
    }
}

function dragElement(element) {
/* Nested function to track selected element, mouse dragging, and mouse release. If any elements
are stored in element.src_lines or element.dst_lines, their top/left properties are moved the
same way as the selected element.
Parameters
----------
element : HTMLElement
    Element to which to apply dragging functionality.
*/
    // To keep track of previous position and distance moved by cursor.
    var oldX = 0, oldY = 0;
    var moveX = 0, moveY = 0;
    element.onmousedown = grab;
    element.class = NODECLASS;
    element.ondblclick = openPanel;
    element.move = move;
    element.del = del;
    lastSelected = element;  // last created == last selected
    function grab(e) {
        // Selects the element and adds the onmousemove function to have the element follow the cursor.
        e = e || window.event();
        e.preventDefault();
        lastSelected = element;  // marks the current element as grabbed
        oldX = e.clientX;  // Current mouse position
        oldY = e.clientY;  // Current mouse position
        document.onmouseup = release;  // while mouse is held, keep dragging
        document.onmousemove = drag;
    }
    function drag(e) {
        // Updates position of the currently-grabbed element
        e = e || window.event;
        e.preventDefault();
        // Get move vector, update current position
        moveX = e.clientX - oldX;
        moveY = e.clientY - oldY;
        oldX = e.clientX;
        oldY = e.clientY;
        move(moveX, moveY);
    }

    function move(moveX, moveY){
        // set new pos
        element.style.top = (element.offsetTop + moveY) + "px";
        element.style.left = (element.offsetLeft + moveX) + "px";
        // move attached lines, if any
        for (let idx = 0; idx < element.src_lines.length; idx++){
          el_to_move = document.getElementById(element.src_lines[idx]);
          el_to_move.srcX += moveX;
          el_to_move.srcY += moveY;
          el_to_move.update_pos();
        }
        for (let idx = 0; idx < element.dst_lines.length; idx++){
          el_to_move = document.getElementById(element.dst_lines[idx]);
          el_to_move.dstX += moveX;
          el_to_move.dstY += moveY;
          el_to_move.update_pos();
        }
        // move label, if any
        textDiv = document.getElementById(element.label);
        textDiv.style.top = addToPx(textDiv.style.top, moveY);
        textDiv.style.left = addToPx(textDiv.style.left, moveX);
    }
    function del(){
        // delete this element and associated elemeents
        for(let idx=element.dst_lines.length-1; idx>=0; idx--){
          document.getElementById(element.dst_lines[idx]).del();
        }
        for(let idx=element.src_lines.length-1; idx>=0; idx--){
          document.getElementById(element.src_lines[idx]).del();
        }
        // delete text label
        textDiv = document.getElementById(element.label);
        textDiv.remove();
        // remove self
        idx = nodeIds.indexOf(element.id);
        nodeIds.splice(idx, 1);
        element.remove();
        return;

    }


    function release() {
        // releases grabbed element; will no longer be grabbed
        document.onmouseup = null;
        document.onmousemove = null;
   }

   function openPanel(e){
     // Close other panel if it's already opened
     if (lastOpenPanel != null){
       lastOpenPanel.remove();
     }
     // Don't open a panel if node has no parameters
     if (element.params.length == 0){
       return;
     }
     lastSelected = null;  // prevent character deletion from
     const panel = document.createElement("div");
     panel.class = PANELCLASS;
     panel.id = PANELCLASS + SEPSTR + getNodeNumFromId(element.id);
     panel.nodeId = element.id;
     applyClassCSS(panel);
     // Create table for user input
     tbl = document.createElement('table');
     tbl.style.width = '100px';
     tbl.style.border = '5px solid black';
     // Name for display / input selection later
     tr = tbl.insertRow();
     td = tr.insertCell();
     td.appendChild(document.createTextNode('Node label'));
     td = tr.insertCell();
     inp = document.createElement('input');
     inp.defaultValue = element.labelText;
     inp.id = NAMECLASS + SEPSTR + getNodeNumFromId(element.id);
     td.appendChild(inp);
     for (var idx=0; idx<element.params.length; idx++) {
       // Get labels / default values, if any
       parName = element.params[idx];
       parVal = element.paramsValues[parName];
       parType = element.paramsTypes[parName];  // todo: use for validation
       // Create row; one per parameter
       tr = tbl.insertRow();
       td = tr.insertCell();
       td.appendChild(document.createTextNode(parName));
       td = tr.insertCell();
       inp = document.createElement('INPUT');
       if (parVal != null){
         inp.defaultValue = parVal;
       }
       inp.id = INPUTCLASS + SEPSTR + parName;
       inp.parType = parType;
       td.appendChild(inp);
     }

     // Create save/close buttons
     tr = tbl.insertRow();
     td = tr.insertCell();
     btn = document.createElement("button");
     btn.onclick = closePanel;
     btn.panelId = panel.id;
     btn.innerHTML = "&times; Cancel"
     td.appendChild(btn);
     td = tr.insertCell();
     btn = document.createElement("button");
     btn.onclick = saveParams;
     btn.nodeId = element.id;
     btn.panelId = panel.id;
     btn.innerHTML = "&#10003; Save"
     td.appendChild(btn);

     panel.appendChild(tbl);
     document.body.appendChild(panel);
     lastOpenPanel = panel;
   }
}

function closePanel(event) {
  btn = event.target;
  panel = document.getElementById(btn.panelId);
  panel.remove();
  return
}
function saveParams(event) {
  // Save parameters that are in the panel and store them in the node.
  btn = event.target;
  panelId = btn.panelId;
  nodeId = btn.nodeId;
  panel = document.getElementById(panelId);
  node = document.getElementById(panel.nodeId);
  tbl = panel.children[0];
  var params = new Object();
  for (var i=0, row; row = tbl.rows[i]; i++){
    for (var j=0, cell; cell = row.cells[j]; j++){
      if (cell.children.length > 0 && cell.children[0].id.includes(INPUTCLASS)){  // only select boxes of this class
        txtBox = cell.children[0];
        param = txtBox.id.substring(INPUTCLASS.length + 1, txtBox.id.length);
        if (txtBox.value.length > 0){
          params[param] = txtBox.value;
        }
        else if (txtBox.value.length == 0){
          params[param] = null;
        }
      }
      else if(cell.children.length > 0 && cell.children[0].id.includes(NAMECLASS)){
        // display name for the node
        txtBox = cell.children[0];
        name = txtBox.value;
        updateLabel(node, name);
      }
    }
  }
  node.paramsValues = params;
  closePanel(event);
  return
}

function createDraggable(img="static/test2.png", name, module, params, analysisNodeType="simo", id) {
/* Creates a draggable image with an imagemap and adds them to the DOM.


*/
  if(id == null){
    sample_id = nodeCount++;
  }
  else{
    sample_id = id;
    nodeCount = parseInt(id) + 1;
  }
  // create elements
  draggable = createDraggableImage(sample_id, HORIZONTAL_POSITION_INITIAL, VERTICAL_POSITION_INITIAL, img, name);
  draggable.analysisNodeType = analysisNodeType;
  imagemap = createDraggableImageMap(sample_id);

  document.body.appendChild(draggable);
  document.body.appendChild(imagemap);
  var param_names = new Array();
  draggable.paramsValues = new Object();
  draggable.paramsTypes = new Object();
  params = parseParams(params);
  for (const p in params){
    param_names.push(p);
    draggable.paramsTypes[p] = params[p][0];
    draggable.paramsValues[p] = params[p][1];
  }
  draggable.title = param_names;
  draggable.params = param_names;

  draggable.module = module;
  draggable.procName = name;
  nodeIds.push(draggable.id);
  return draggable;
}

function createDraggableImage(id, x = 200, y = 100, img="static/test2.png", name="node") {
/* Returns a draggable image element at the specified coordinates using the input image.
Parameters
----------
id : int
    Integer to use for the id; id is of the form [DRAGSTR][SEPSTR][id]. This id should also
    be used with createDraggableImageMap to match the two.
x, y : int
    Position on the page where to add the image.
img : String
    Path to the image to use for the image
Returns
-------
element
    Handle of the created draggable element.
*/
  const draggable = document.createElement("img");
  draggable.class = NODECLASS;
  applyClassCSS(draggable);

  draggable.src = img;
  draggable.src_img = img;
  draggable.deg = 0;
  draggable.draggable = "true";
  draggable.width = "81";
  draggable.height = "81";
  draggable.id = DRAGSTR + SEPSTR + id;
  draggable.useMap = '#' + IMAPSTR + SEPSTR + id;

  draggable.style.left = x.toString() + 'px';
  draggable.style.top = y.toString() + 'px';

  const textDiv = document.createElement("div");
  textDiv.class = TEXTCLASS;
  textDiv.innerHTML = name;
  applyClassCSS(textDiv);
  textDiv.style.left = addToPx(draggable.style.left, 0);
  textDiv.style.top = addToPx(draggable.style.top, -10);
  textDiv.style.zIndex = '2';
  textDiv.id = TEXTCLASS + SEPSTR + id;
  document.body.appendChild(textDiv);
  draggable.label = textDiv.id;
  draggable.labelText = name;

  draggable.src_lines = new Array();  // edges that start at this node
  draggable.dst_lines = new Array();  // edges that end at this node
  dragElement(draggable);
  return draggable;
}

function updateLabel(node, newname){
  textDiv = document.getElementById(node.label);
  node.labelText = newname;
  textDiv.innerHTML = node.labelText;

}

function createDraggableImageMap(id) {
/* Creates an ImageMap to be overlaid on a draggable image with a corresponding ID.
Parameters
----------
id : int
    ID used to create the draggable image.
Returns
-------
element
    Handle for the element.
*/
  image_map = document.createElement("map");
  image_map.name = IMAPSTR + SEPSTR + id;
  image_map.id = image_map.name;

  // Input area specs
  image_area = document.createElement("area");
  image_area.id = AREASTR + SEPSTR + id + SEPSTR + 'in' + SEPSTR + '0';
  image_area.shape = "rect";
  image_area.coords = "30,0,50,20";
//  image_area.onmouseup = function() {alert('thing!')};
  image_area.onmouseenter = hoverDot;
  image_area.onmouseout = unhoverDot;
  image_map.appendChild(image_area);

  // Output area specs
  image_other_area = document.createElement("area");
  image_other_area.id = AREASTR + SEPSTR + id + SEPSTR + 'out' + SEPSTR + '0';
  image_other_area.shape = "rect";
  image_other_area.coords = "30,50,50,70";
  image_other_area.onmousedown = lineCreator;
  image_other_area.onmouseenter = hoverDot;
  image_other_area.onmouseout = unhoverDot;
  image_map.appendChild(image_other_area);
  return image_map;
}

function lineCreator(event) {
/* Handler to create a scalable line centered on the target of the event, which is intended to
 be an area of an imagemap on a draggable image.
Parameters
----------
event
    Event to handle.
*/
    // get target (imagemap area)
    area = event.target;
    coords = parseIntList(area.coords);
    area_width_center = (coords[2] - coords[0])/2;
    // get id of map; these are made to correspond to draggable node
    mapid = area.parentElement.id;
    node_num = mapid.split(SEPSTR).pop();
    draggable_id = DRAGSTR + SEPSTR + node_num;
    draggable = document.getElementById(draggable_id);
    dragLeft = parseInt(draggable.style.left.split('px')[0]);
    dragTop = parseInt(draggable.style.top.split('px')[0]);
    srcX = dragLeft + parseInt(coords[0]) + area_width_center;  // center of the area
    srcY = dragTop + parseInt(coords[3]);  // bottom of the area

    // Get line id; need to make sure that it is unique
    line_count = draggable.src_lines.length;
    line_num = node_num + SEPSTR + line_count.toString();  // id for this line
    line_id = LINESTR + SEPSTR + line_num;
    while (draggable.src_lines.includes(line_id)){
      line_count += 1;
      line_num = node_num + SEPSTR + line_count.toString();  // id for this line
      line_id = LINESTR + SEPSTR + line_num;

    }
//    line.srcNode = draggable_id;
    line = createLine(line_num, srcX, srcY, event.clientX, event.clientY, draggable_id);
    line_id = line.id;

//    draggable.src_lines.push(line_id);
}
function createLine(line_num, srcX=0, srcY=0, dstX=10, dstY=100, srcNodeId="") {
/* Creates a div element that acts as a connector between nodes. The div can be dragged and follows
the cursor until released. If released over the input of a node, the line attaches itself to the node.
*/
    const line = document.createElement("div");
//    line.id = LINESTR + SEPSTR + line_num.toString();
    line.class = "connector";
    applyClassCSS(line);

    line.onmousedown = function() {lastSelected = line};
    line.srcNode = srcNodeId;
    line.dstNode = NaN;
    // Set position params
    line.srcX = srcX;
    line.srcY = srcY;
    line.dstX = dstX;
    line.dstY = dstY;

    // Compute position parameters
    height = Math.sqrt((dstX-srcX)**2 + (dstY-srcY)**2);
    rot_angle = getRadAngleFromVertical(dstY-srcY, dstX-srcX);
    rot_deg = radToDeg(rot_angle);
    line.deg = rot_deg;

    // Apply params to div style
    line.style.height = height.toString() + "px";
    line.style.left = (line.srcX).toString() + "px";
    line.style.top = (line.srcY).toString() + "px";
    line.style.transform = "rotate(" + line.deg.toString() + "deg)";

    // Apply movement properties
    stretchLine(line);
    document.body.appendChild(line);
    return line;
}

function stretchLine(element) {
/* Nested function to stretch a line div to follow the cursor or nodes it is attached to.
Parameters
----------
element
    Element handle for the line div.
*/
    // Position parameters
    var moveX = 0, moveY = 0;
    element.deg = 0;
    lastSelected = element;
    document.onmousemove = drag;
    document.onmouseup = release;
    element.update_pos = update_pos;
    element.del = del;
    lastSelected = element;  // last created == last selected

    function drag(e) {
        // updates position of the currently-grabbed element
        e = e || window.event;
        e.preventDefault();
        element.dstX = e.clientX;
        element.dstY = e.clientY;
        update_pos(e);
    }

//    function update_pos(event) {
//      event = event || window.event;
//      event.preventDefault();
    function update_pos() {
      // Position parameters determined by src
      element.style.left = (element.srcX).toString() + "px";
      element.style.top = (element.srcY).toString() + "px";

      // Height determined by src, dist
      height = Math.sqrt((element.dstX - element.srcX)**2 + (element.dstY - element.srcY)**2);
      height = Math.round(height);
      element.height = height.toString() + "px";
      element.style.height = height.toString() + "px";

      // Rotation determiend by src, dst; stored in .deg
      deltaX = element.dstX - element.srcX;
      deltaY = element.dstY - element.srcY;
      rad_angle = getRadAngleFromVertical(deltaY, deltaX);
      rot_deg = radToDeg(rad_angle);
      element.deg = rot_deg;
      element.style.transform = "rotate(" + element.deg.toString() + "deg)";
    }

    function del() {
        // Go to srcNode, remove line id from src list
        srcNode = document.getElementById(element.srcNode);
        newSrcLines = new Array();
        newDstLines = new Array();
        for (let idx=0; idx<srcNode.src_lines.length; idx++) {
          if (srcNode.src_lines[idx] != element.id){
            newSrcLines.push(srcNode.src_lines[idx]);
          }

        }
        srcNode.src_lines = newSrcLines;

        // Go to dstNode, remove line id from dst list
        dstNode = document.getElementById(element.dstNode);
        for (let idx=0; idx<dstNode.dst_lines.length; idx++) {
          if (dstNode.dst_lines[idx] != element.id){
            newDstLines.push(dstNode.dst_lines[idx]);
//            dstNode.dst_lines = dstNode.dst_lines.splice(idx, idx);
//            break;
          }
        }
        dstNode.dst_lines = newDstLines;
        element.remove();
    }


    function release(e) {
        // releases grabbed element; will no longer be grabbed
        // find nearest draggable; see if we can attach
//        node_id = element.id.split(SEPSTR)[1];
        node_id = element.srcNode.split(SEPSTR)[1];
        drag_list = document.elementsFromPoint(e.clientX, e.clientY);
        var foundMatch = 0;
        for (let idx=0; idx < drag_list.length; idx++){
          drag_el = drag_list[idx];
          if(drag_el.id.startsWith(DRAGSTR)){
            second_node_id = drag_el.id.split(SEPSTR)[1];
            if(second_node_id == node_id){
              // skip if trying to attach to self
              continue;
            }
            // get draggable coords; get area 0 coords
            imap_id = IMAPSTR + SEPSTR + second_node_id;
            imap = document.getElementById(imap_id);

            // for each area in the imagemap, check if it is an input
            // compute the distance between the current position and the center of the area
            // connect to closest
            closest = Infinity;
            closest_area = NaN;
            closest_center = NaN;
            for (let area_idx=0; area_idx<imap.childNodes.length; area_idx++){
              area = imap.childNodes[area_idx]
              if (area.id.includes(SEPSTR + 'in' + SEPSTR)){
                // is an input
                area_coords = area.coords;
                p_center = getCenterOfRect(area_coords);
                dist = computeL2Dist(p_center, [e.clientX, e.clientY]);
                if (dist < closest){
                    closest = dist;
                    closest_area = area;
                    closest_center = p_center;
                }
              }
            }
            if (!isNaN(closest_area)) {
              // connect to center
              element.dstNode = drag_el.id;
              srcNode = document.getElementById(element.srcNode);

              // create line id
              srcNum = getNodeNumFromId(srcNode.id);
              dstNum = getNodeNumFromId(drag_el.id);
              lineId = LINESTR + SEPSTR + srcNum + SEPSTR + dstNum;
              // check if this id exists already
              otherLine = document.getElementById(lineId);
              // todo: fix this
              if (otherLine != null){  // a line connecting the two nodes already exists; don't connect
                break;
              }
              element.id = lineId;
              element.dstX = p_center[0] + drag_el.offsetLeft;
              element.dstY = p_center[1] + drag_el.offsetTop;

              srcNode = document.getElementById(element.srcNode);
              srcNode.src_lines.push(element.id);
              drag_el.dst_lines.push(element.id);
              element.update_pos(e);
              foundMatch = 1;
            }
            break;
          }
        }
        if (foundMatch == 0) {
          element.remove();
        }
        document.onmouseup = null;
        document.onmousemove = null;
   }
}

function getRadAngleFromVertical(deltaY, deltaX){
  // Gets the angle of rotation for a displacement relative to the vertical axis.
  atan = Math.atan(deltaY / deltaX)
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

function addToPx(px_str, value) {
  // extract value; add value
  px_val = parseInt(px_str.split('px')[0]) + value;
  return px_val.toString() + "px";
}

function parseIntList(s, sep=',') {
  slist = s.split(sep);
  arr = new Array(slist.length);
  for (let idx=0; idx<arr.length; idx++){
    arr[idx] = parseInt(slist[idx]);
  }
  return arr;
}

function getCenterOfRect(coords){
// finds the center of the rectangle specified by coords
// expected format: x0,y0,x1,y1
  pos = parseIntList(coords);
  cx = (pos[0] + pos[2])/2;
  cy = (pos[1] + pos[3])/2;
  return [cx, cy];
}

function computeL2Dist(p0, p1) {
// computes the L2 distance between points p0 and p1, where each is an array of the form [x, y]
return Math.sqrt((p1[0]-p0[0])**2 + (p1[1]-p0[1])**2);
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
var class_to_check = element.class;
if (cls.length != 0) {
  class_to_check = cls;
}
// Find the matching rule
for (const sheet of document.styleSheets) {
    for (const r of sheet.rules){
        splitText = r.selectorText.replace(" ","").split(",");
        if (splitText.includes("." + class_to_check)){
            element.style = r.style.cssText;
        }
    }
  }
}

function hoverDot(event) {
// Makes a dot appear near the center of the target. This is intended to highlight node connection points.
  area = event.target;
  area_coords = area.coords;
  coords = parseIntList(area.coords);
  p_center = getCenterOfRect(area_coords);

  // get the dot element, or create if it doesn't exist
  el = document.getElementById("bluedot");
  if (el == null){
    el = document.createElement("span");
    el.id = "bluedot";
    el.class = "hoverdot";
    applyClassCSS(el);
    document.body.appendChild(el);
  }
  // distinguish between input and output points
  var vertOffset = 0
  if (area.id.includes("in")){
    el.style.backgroundColor = "#0000bb";
    vertOffset = coords[1] - 3 ; // move the dot slightly higher for visibility; input assumed to be above
  }
  else {
    el.style.backgroundColor = "#bb0000";
    vertOffset = p_center[1]; // move dot to center of area.
  }
  el_height = parseInt(el.style.height.split("px")[0]);
  // move dot to appropriate location
  draggable = getDraggableFromImapID(area.parentElement.id);
  el.style.top = (draggable.offsetTop + vertOffset).toString() + "px";
  el.style.left = (draggable.offsetLeft + p_center[0]-el_height/2).toString() + "px";
  el.style.opacity = 0.5;
  el.style.zIndex = 1;
}

function unhoverDot(event) {
// hide dot
  el = document.getElementById("bluedot");
  if (el == null){
    return;
  }
  el.style.opacity = 0;
  el.style.zIndex = -1;

}

function getDraggableFromImapID(id) {
  node_num = id.split(SEPSTR).pop();
  draggable_id = DRAGSTR + SEPSTR + node_num;
  return document.getElementById(draggable_id);
}

function parseParams(s) {
// This function parses a string-ified Python list
  // remove brackets
  if(typeof s == "string"){
    s = s.replaceAll("'", '"');
    s = s.replaceAll("None", null);
    params = JSON.parse(s);
    return params;
  }
  else{
    return s;
  }

}

function getNodeNumFromId(id){
  return id.split(SEPSTR).pop();
}

function exportPipeline(pipelineName, pipelineDestination) {
// This function pulls together all nodes, their properties, and connections
  isDestProject = pipelineDestination.startsWith('project:');
  dest = pipelineDestination.substring(pipelineDestination.indexOf(':')+1);

  var node_export = new Array();
  for(var i=0; i<nodeIds.length; i++){
    el = document.getElementById(nodeIds[i]);
    el.nodeId = el.id;
    node_export.push(el);
  }
  var pipelineSummary = new Object();
  pipelineSummary['name'] = pipelineName;
  pipelineSummary['nodes'] = node_export;
  pipelineSummary['isDestProject'] = isDestProject;
  pipelineSummary['saveId'] = dest;
    fetch("/pipeline", {
        method: "POST",
        headers: {
            "Accept": "application/json",
            "Content-Type": "application/json"
        },
        body: JSON.stringify(pipelineSummary)
    });
}

function openSavePanel(projDestsJSON, userDestsJSON) {
  projDestsJSON = projDestsJSON.replaceAll("'", '"');
  projDestsJSON = projDestsJSON.replaceAll("None", null);
  projDests = JSON.parse(projDestsJSON);

  userDestsJSON = userDestsJSON.replaceAll("'", '"');
  userDestsJSON = userDestsJSON.replaceAll("None", null);
  userDests = JSON.parse(userDestsJSON);

  const panel = document.createElement("div");
  panel.class = PANELCLASS;
  panel.id = PANELCLASS + SEPSTR + "save";
  applyClassCSS(panel);

  tbl = document.createElement('table');
  tbl.style.width = '100px';
  tbl.style.border = '5px solid black';

  // Need: pipeline name, associated project
  tr = tbl.insertRow();
  td = tr.insertCell();
  inpName = document.createElement("INPUT");
  td.appendChild(inpName);

  td = tr.insertCell();
  drop = document.createElement("select");
  drop.style.maxWidth = '100px';
  option = document.createElement("option");
  option.innerHTML = "Destination";
  option.disabled = true;
  option.selected = true;
  drop.appendChild(option);

  option = document.createElement("option");
  option.innerHTML = "Userspace: " + userDests['user'];
  option.value = 'user:' + userDests['id'];
  drop.appendChild(option);
  for (var proj in projDests){
    option = document.createElement("option");
    var projDisplay = proj;
    if (proj.length > 20){
      projDisplay = proj.substring(0,17) + "...";
    }
    option.innerHTML = "Project: " + projDisplay;
    option.value = 'project:' + projDests[proj];
    drop.appendChild(option);
  }
  td.appendChild(drop);

  td = tr.insertCell();
  btn = document.createElement("button");
  btn.innerHTML = "Save";
  btn.onclick = function() {exportPipeline(inpName.value.toString(), drop.value.toString() ); panel.remove();}
  td.appendChild(btn);
  panel.appendChild(tbl);

  document.body.appendChild(panel);
}

async function getAnalysis(){
  analysis_selection = document.getElementById("analyses");
  id_analysis = analysis_selection.value;
  var analysis_spec = new Object();
  analysis_spec["loadAnalysis"] = id_analysis;

  let response = await fetch(window.location.href, {
                          method: "POST",
                          headers: {
                                "Accept": "application/json",
                                "Content-Type": "application/json"},
                          body: JSON.stringify(analysis_spec)});
  let data = await response.json();
  nodes = data["nodes"];
  clearNodes();
  var nodeObject = new Object();
  var nodesWithoutDepth = new Array();
  for (let idx=0; idx<nodes.length; idx++){
    node = nodes[idx];
    var load_params = new Object();
    for (let param_idx=0; param_idx<node["params"].length; param_idx++){
      paramName = node["params"][param_idx];
      paramType = node["paramsTypes"][paramName];
      paramValue = node["paramsValues"][paramName];
      load_params[paramName] = [paramType, paramValue];
    }
    id = getNodeNumFromId(node["nodeId"]);
    loadedNode = createDraggable(nodes[idx]["src_img"], nodes[idx]["labelText"], nodes[idx]["module"], load_params, node["analysisNodeType"], id);
    if(node["analysisNodeType"] == "input"){
      loadedNode.depth = 0;
    }
    else{
      loadedNode.depth = -1;
      nodesWithoutDepth.push(loadedNode);
    }
    loadedNode.id = node["nodeId"];

//    Object.assign(loadedNode, node["nodeId"]);
    loadedNode.dst_lines = node["dst_lines"];
    loadedNode.src_lines = node["src_lines"];
    nodeNum = getNodeNumFromId(node["nodeId"]);
    nodeObject[nodeNum] = loadedNode;
  }

  // make connectors
  for(var nodeIdx in nodeObject){
    currentNode = nodeObject[nodeIdx];
    // go through src
    currentNodeNum = getNodeNumFromId(currentNode.id);
    srcAreaPos = getAreaCenter(currentNode, idx=0, areaType="out");
    srcNodeLeft = parseInt(currentNode.style.left.split("px")[0]);
    srcNodeTop = parseInt(currentNode.style.top.split("px")[0]);
    for(var srcIdx in currentNode.src_lines){
        srcString = currentNode.src_lines[srcIdx];
        nextNodeNum = getNodeNumFromSrc(srcString);
        nextNode = getNodeFromNum(nextNodeNum);

        // We want to connect the output of currentNode to the input of nextNode
        nextAreaPos = getAreaCenter(nextNode, idx=0, areaType="in");
        nextNodeLeft = parseInt(nextNode.style.left.split("px")[0]);
        nextNodeTop = parseInt(nextNode.style.top.split("px")[0]);
        line = createLine("", srcAreaPos[0]+srcNodeLeft, srcAreaPos[1]+srcNodeTop, nextAreaPos[0]+nextNodeLeft, nextAreaPos[1]+nextNodeTop, currentNode.id);
        line.id = srcString;
        line.dstNode = nextNode.id;
        document.onmousemove = null;
        document.onmouseup = null;
    }

  }
  // Below this is aesthetics.
  // determine depth of each node
  while(nodesWithoutDepth.length > 0){
    var spliceIdxList = new Array();
    for(let spliceIdx = 0; spliceIdx < nodesWithoutDepth.length; spliceIdx++){
      node = nodesWithoutDepth[spliceIdx];
      // examine previous nodes; if they all have depth, assign max depth+1 to this node
      var maxDepth = -1;
      var doBreak = false;
      for(let previousNodeIdx = 0; previousNodeIdx < node["dst_lines"].length; previousNodeIdx++){
        // get num of previous node
        var prevNum = getNodeNumFromDst(node["dst_lines"][previousNodeIdx]);
        var prevNode = nodeObject[prevNum];
        if(prevNode.depth == -1){
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
  var xposByDepth = new Object();
  for(let idx=0; idx<maxDepth+2; idx++){
    xposByDepth[idx] = 0;
  }
  for(var n in nodeObject){
    var node = nodeObject[n];
    node.xpos = xposByDepth[node.depth];
    xposByDepth[node.depth]+=1;
  }

  // move nodes to position
  for(idx in nodeObject){
    node = nodeObject[idx];
    xpos = node.xpos * HORIZONTAL_SPACING + Math.random() * HORIZONTAL_JIGGLE;
    ypos = node.depth * VERTICAL_SPACING + Math.random() * VERTICAL_JIGGLE;
    node.move(xpos, node.depth * VERTICAL_SPACING + Math.random() * VERTICAL_JIGGLE);
  }

}

function getNodeNumFromDst(dstString){
/* Returns the node number of the node leading to the current node, as determined by the dstString.
Expected structure: LINESTR + SEPSTR + PREVIOUS + SEPSTR + CURRENT
*/
return dstString.split(SEPSTR)[1]
}

function getNodeNumFromSrc(srcString){
/* Returns the node number of the node downstream from the current node, as determined by the srcString.
Expected structure: LINESTR + SEPSTR + CURRENT + SEPSTR + NEXT
*/
return srcString.split(SEPSTR)[2];
}

function getNodeFromNum(num){
/* Returns the node specified by the number */
return document.getElementById(DRAGSTR + SEPSTR + String(num));
}

function getImageMapFromNode(node, areaIdx=0, areaType="in"){
/* Gets the imagemap for the node; if there are multiple maps, return the one specified by inputIdx
node : object representing pipeline nodes
areaIdx : index for which map to return
areaType : "in" or "out"; specify the type of area being identified
*/
  nodeNum = getNodeNumFromId(node.id);
  areaMapId = AREASTR + SEPSTR + nodeNum + SEPSTR + areaType + SEPSTR + String(areaIdx);
  return document.getElementById(areaMapId);
}

function getAreaCenter(node, idx=0, areaType="in"){
  // gets the X, Y coordinates for the center of an area
  area = getImageMapFromNode(node, idx, areaType);
  coords = parseIntList(area.coords);
  centerX = Math.floor((coords[0] + coords[2])/2);
  centerY = Math.floor((coords[1] + coords[3])/2);
  return [centerX, centerY]
}

function clearNodes(){
  // removes all nodes from the page
  for(let idx=nodeIds.length-1; idx>=0; idx--){
    nodeId = nodeIds[idx];
    node = document.getElementById(nodeId);
    node.del();
  }
  nodeCount = 0;
}