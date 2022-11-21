document.addEventListener("keydown", delTest, false);
var lastSelected = null;
var nodeCount = 0;
const SEPSTR = '-';
const DRAGSTR = "draggable";
const LINESTR = "line";
const IMAPSTR = "imagemap";
const AREASTR = "area";

function delTest(e) {
    if (event.keyCode == 46) {
        lastSelected.remove();
    }
}

function dragElement(element) {
    var moveX = 0, moveY = 0;
    var oldX = 0, oldY = 0;
    element.onmousedown = grab;
    lastSelected = element;  // last created == last selected
    function grab(e) {
        // marks the current element as grabbed
        e = e || window.event();
        e.preventDefault();
        lastSelected = element;
        oldX = e.clientX;
        oldY = e.clientY;
        document.onmouseup = release;  // while mouse is held, keep dragging
        document.onmousemove = drag;
    }
    function drag(e) {
        // updates position of the currently-grabbed element
        e = e || window.event;
        e.preventDefault();
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
          el_to_move.style.top = addToPx(el_to_move.style.top, moveY);
          el_to_move.style.left = addToPx(el_to_move.style.left, moveX);
        }
    }
    function release() {
        // releases grabbed element; will no longer be grabbed
        document.onmouseup = null;
        document.onmousemove = null;
   }
}

function createDraggable(x = 0, y = 0, img = "static/test2.png") {
  sample_id = nodeCount++;
  // create elements
  draggable = createDraggableImage(sample_id);
  imagemap = createDraggableImageMap(sample_id);

//  document.body.appendChild(draggable);
  document.body.appendChild(imagemap);
  dragElement(draggable);
}

function createDraggableImage(id, x = 0, y = 0, img="static/test2.png") {
  const draggable = document.createElement("img");
  draggable.class = "dragged";
  draggable.src = img;
  draggable.deg = 0;
  draggable.draggable = "true";
  draggable.style = "position:absolute;cursor:move;";
  draggable.width = "81";
  draggable.height = "81";
  draggable.id = DRAGSTR + SEPSTR + id;
  draggable.useMap = '#' + IMAPSTR + SEPSTR + id;
  document.body.appendChild(draggable);
  draggable.style.top = "150px";
  draggable.style.left = "50px";
  draggable.src_lines = new Array();  // edges that start at this node
  draggable.dst_lines = new Array();  // edges that end at this node
  return draggable;
}

function createDraggableImageMap(id) {
  image_map = document.createElement("map");
  image_map.name = IMAPSTR + SEPSTR + id;
  image_map.id = image_map.name;
  image_area = document.createElement("area");
  image_area.id = AREASTR + SEPSTR + id + SEPSTR + 'in' + SEPSTR + '0';
  image_area.shape = "rect";
  image_area.coords = "30,0,50,20";
  image_area.onmouseup = function() {alert('thing!')};
  image_map.appendChild(image_area);

  image_other_area = document.createElement("area");
  image_other_area.id = AREASTR + SEPSTR + id + SEPSTR + 'out' + SEPSTR + '1';
  image_other_area.shape = "rect";
  image_other_area.coords = "30,50,50,70";
  image_other_area.onmousedown = lineCreator;
  image_map.appendChild(image_other_area);
  return image_map;
}

function lineCreator(event) {
    // get imagemap area
    area = event.target;
    coords = area.coords.split(",");
    area_width = (parseInt(coords[2]) - parseInt(coords[0]))/2;
    // get id of map; these are made to correspond to draggable node
    mapid = area.parentElement.id;
    node_num = mapid.split(SEPSTR).pop();
    draggable_id = DRAGSTR + SEPSTR + node_num;
    draggable = document.getElementById(draggable_id);
    dragLeft = parseInt(draggable.style.left.split('px')[0]);
    dragTop = parseInt(draggable.style.top.split('px')[0]);
    srcX = dragLeft + parseInt(coords[0]) + area_width;
    srcY = dragTop + parseInt(coords[3]);
    line_num = node_num + SEPSTR + (draggable.src_lines.length).toString();
    line_id = createLine(line_num, srcX, srcY, event.clientX, event.clientY);
    draggable.src_lines.push(line_id);
}
function createLine(line_num, srcX=0, srcY=0, dstX=10, dstY=100) {
    const line = document.createElement("div");
    line.id = LINESTR + SEPSTR + line_num.toString();
    line.class = "line";
    line.srcX = srcX;
    line.srcY = srcY;
    line.dstX = dstX;
    line.dstY = dstY;
    line.style.backgroundColor = 'black';
    height = Math.sqrt((dstX-srcX)**2 + (dstY-srcY)**2);
    rot_angle = getRadAngleFromVertical(dstY-srcY, dstX-srcX);
    rot_deg = radToDeg(rot_angle);
    line.deg = rot_deg;
    line.style.height = height.toString() + "px";
    line.style.width = "5px";

    line.style.left = (line.srcX).toString() + "px";
    line.style.top = (line.srcY).toString() + "px";
    line.style.position = "absolute";
    line.style.transformOrigin = "50% 0%"
    line.style.transform = "rotate(" + line.deg.toString() + "deg)";

    el = document.getElementById(DRAGSTR + SEPSTR + node_num.toString());
    document.body.appendChild(line);
    stretchLine(line);
    return line.id;
}

function stretchLine(element) {
    var moveX = 0, moveY = 0;
    element.deg = 0;
    lastSelected = element;
    document.onmousemove = drag;
    document.onmouseup = release;
    lastSelected = element;  // last created == last selected
    function grab(e) {
        // marks the current element as grabbed
        e = e || window.event();
        e.preventDefault();
        lastSelected = element;
        document.onmouseup = release;  // while mouse is held, keep dragging
        document.onmousemove = drag;
    }
    function drag(e) {
        // updates position of the currently-grabbed element
        e = e || window.event;
        e.preventDefault();
        element.dstX = e.clientX;
        element.dstY = e.clientY;
        moveX = element.dstX - element.srcX;
        moveY = element.dstY - element.srcY;
        rad_angle = getRadAngleFromVertical(moveY, moveX);
        rot_deg = radToDeg(rad_angle);
        element.deg = rot_deg;
        height = Math.sqrt((element.dstX - element.srcX)**2 + (element.dstY - element.srcY)**2);
        element.height = Math.round(height).toString() + "px";
        element.style.height = Math.round(height).toString() + "px";
        element.style.transform = "rotate(" + element.deg.toString() + "deg)";
//        element.style.scale = ";scale(" + (height/20).toString() + ", 1)";
    }
    function release(e) {
        // releases grabbed element; will no longer be grabbed
        // find nearest draggable; see if we can attach
        node_id = element.id.split(SEPSTR)[1];
        drag_list = document.elementsFromPoint(e.clientX, e.clientY);
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
            console.log(imap);

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
              element.dstX = p_center[0];
              element.dstY = p_center[1];
            }


            break;
          }
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

function computeL2Dist(p0, p1){
// computes the L2 distance between points p0 and p1, where each is an array of the form [x, y]
return Math.sqrt((p1[0]-p0[0])**2 + (p1[1]-p0[0])**2);
}