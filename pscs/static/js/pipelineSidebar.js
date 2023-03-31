var sidebarIsOpen = false;

function toggleSidebar() {
  sidebar = document.getElementById("sidebar");
  if(!sidebarIsOpen){
    // sidebar is closed; open it
    applyClassCSS(sidebar, "sidebar-open");
    // move toggle button
    toggleBtn = document.getElementById("sideToggle");
    toggleBtn.style.left = addTwoPx(getComputedStyle(toggleBtn).left, getComputedStyle(sidebar).width);
    toggleBtn.innerHTML = "&larr;";
    sidebarIsOpen = true;
  }
  else {
    applyClassCSS(sidebar, "sidebar-closed");
    sidebarIsOpen = false;
    toggleBtn = document.getElementById("sideToggle");
//    toggleBtn.style.left = addTwoPx(getComputedStyle(toggleBtn).left, getComputedStyle(sidebar).width);
    toggleBtn.style.left = addToPx(getComputedStyle(toggleBtn).left, -getPxValue(getComputedStyle(sidebar).width));
    toggleBtn.innerHTML = "&rarr;";
  }
  document.body.style.marginLeft = addTwoPx(sidebar.style.left, getComputedStyle(sidebar).width);
}

function getPxValue(pxString){
  return parseInt(pxString.split("px")[0]);
}
function addToPx(pxString, value){
  pxValue = parseInt(pxString.split("px")[0]);
  return (pxValue + value).toString() + "px";
}

function addTwoPx(pxString0, pxString1){
  pxValue0 = parseInt(pxString0.split("px")[0]);
  pxValue1 = parseInt(pxString1.split("px")[0]);
  return (pxValue0 + pxValue1).toString() + "px";
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