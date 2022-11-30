var sidebarIsOpen = false;

function toggleSidebar() {
  sidebar = document.getElementById("sidebar");
  if(!sidebarIsOpen){
    // sidebar is closed; open it
    applyClassCSS(sidebar, "sidebar-open");
    sidebarIsOpen = true;
  }
  else {
    applyClassCSS(sidebar, "sidebar-closed");
    sidebarIsOpen = false;
  }
  document.body.style.marginLeft = sidebar.style.width;
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