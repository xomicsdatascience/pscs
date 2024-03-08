function updateUploaderScripts(scriptArray){
    let s = document.getElementsByClassName("uploader");
    let slength = s.length;
    if(slength !== 0){
        for(let idx=0; idx<slength; idx++){
            document.head.removeChild(s[0]);
        }
    }
    for(let sc of scriptArray){
        addScript(sc, "uploader");
    }
}

function addScript(scriptSrc, scriptClass="uploader"){
    let script = document.createElement("script");
    script.type = "text/javascript";
    script.src = scriptSrc;
    script.classList.add(scriptClass);
    document.head.appendChild(script);
}

function removeScriptsFromDocument(scriptClass){
    let current_scripts = el.getElementsByClassName(scriptClass);
    for (let i = current_scripts.length-1; i >=0; i--){
        if (current_scripts[i].tagName === "SCRIPT"){
            current_scripts[i].parentNode.removeChild(current_scripts[i]);
        }
    }
}