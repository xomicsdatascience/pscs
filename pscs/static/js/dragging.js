// This JS file contains code for dragging, stretching, and reorienting
// HTML elements.

// Global variables
var lastSelected = null;  // currently-selected HTML element
var nodeCount = 0;  // used for creating new node & map IDs

// ID substrings
const SEPSTR = '-';
const DRAGSTR = "draggable";
const LINESTR = "connector";
const IMAPSTR = "imagemap";
const AREASTR = "area";
// Facilitate code reading
const DELKEY = 46;

// For removing elements.
document.addEventListener("keydown", delNode, false);
function delNode(e) {
    if (event.keyCode == DELKEY) {
      len = lastSelected.src_lines.length;
      for (let idx=0; idx<len; idx++){
        remove_id = lastSelected.src_lines.pop();
        line_elem = document.getElementById(remove_id);
        line_elem.del();
      }
      len = lastSelected.dst_lines.length;
      for (let idx=0; idx<len; idx++){
        remove_id = lastSelected.dst_lines.pop();
        line_elem = document.getElementById(remove_id);
        line_elem.del();
      }
      lastSelected.remove();
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
    var currentPos = [0,0];
    var moveVector = [0,0];
    var oldX = 0, oldY = 0;
    var moveX = 0, moveY = 0;
    element.onmousedown = grab;
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
        // set new pos
        element.style.top = (element.offsetTop + moveY) + "px";
        element.style.left = (element.offsetLeft + moveX) + "px";
        // move attached lines, if any
        for (let idx = 0; idx < element.src_lines.length; idx++){
          el_to_move = document.getElementById(element.src_lines[idx]);
          el_to_move.srcX += moveX;
          el_to_move.srcY += moveY;
          el_to_move.update_pos(e);
        }
        for (let idx = 0; idx < element.dst_lines.length; idx++){
          el_to_move = document.getElementById(element.dst_lines[idx]);
          el_to_move.dstX += moveX;
          el_to_move.dstY += moveY;
          el_to_move.update_pos(e);
        }
    }
    function release() {
        // releases grabbed element; will no longer be grabbed
        document.onmouseup = null;
        document.onmousemove = null;
   }
}

function createDraggable(x = 0, y = 0, img = "static/test2.png") {
/* Creates a draggable image with an imagemap and adds them to the DOM.


*/
  sample_id = nodeCount++;
  // create elements
  draggable = createDraggableImage(sample_id);
  imagemap = createDraggableImageMap(sample_id);
  document.body.appendChild(draggable);
  document.body.appendChild(imagemap);

}

function createDraggableImage(id, x = 200, y = 100, img="static/test2.png") {
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
  draggable.class = "dragged";
  applyClassCSS(draggable);

  draggable.src = img;
  draggable.deg = 0;
  draggable.draggable = "true";
  draggable.width = "81";
  draggable.height = "81";
  draggable.id = DRAGSTR + SEPSTR + id;
  draggable.useMap = '#' + IMAPSTR + SEPSTR + id;

  draggable.style.left = x.toString() + 'px';
  draggable.style.top = y.toString() + 'px';

  draggable.src_lines = new Array();  // edges that start at this node
  draggable.dst_lines = new Array();  // edges that end at this node
  dragElement(draggable);
  return draggable;
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
  image_map.appendChild(image_area);

  // Output area specs
  image_other_area = document.createElement("area");
  image_other_area.id = AREASTR + SEPSTR + id + SEPSTR + 'out' + SEPSTR + '1';
  image_other_area.shape = "rect";
  image_other_area.coords = "30,50,50,70";
  image_other_area.onmousedown = lineCreator;
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
    line = createLine(line_num, srcX, srcY, event.clientX, event.clientY);
    line_id = line.id;
    line.srcNode = draggable_id;
//    draggable.src_lines.push(line_id);
}
function createLine(line_num, srcX=0, srcY=0, dstX=10, dstY=100) {
/* Creates a div element that acts as a connector between nodes. The div can be dragged and follows
the cursor until released. If released over the input of a node, the line attaches itself to the node.
*/
    const line = document.createElement("div");
    line.id = LINESTR + SEPSTR + line_num.toString();
    line.class = "connector";
    applyClassCSS(line);
    line.srcNode = NaN;
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

    function update_pos(event) {
      event = event || window.event;
      event.preventDefault();

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
        node_id = element.id.split(SEPSTR)[1];
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
              element.dstX = p_center[0] + drag_el.offsetLeft;
              element.dstY = p_center[1] + drag_el.offsetTop;
              element.dstNode = drag_el.id;
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
return Math.sqrt((p1[0]-p0[0])**2 + (p1[1]-p0[0])**2);
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
const class_to_check = element.class;
if (cls.length != 0) {
  class_to_check = cls;
}
// Find the matching rule
for (const sheet of document.styleSheets) {
    for (const r of sheet.rules){
        if (r.selectorText == "." + class_to_check){
            element.style = r.style.cssText;
        }
    }
  }
}