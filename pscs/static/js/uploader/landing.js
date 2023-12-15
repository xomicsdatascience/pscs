function addScripts(scriptArray){
    let s = document.getElementsByClassName("uploader");
    let slength = s.length;
    if(slength !== 0){
        for(let idx=0; idx<slength; idx++){
            document.head.removeChild(s[0]);
        }
    }
    let script;
    for(let sc of scriptArray){
        script = document.createElement("script");
        script.type = "text/javascript";
        script.src = sc;
        script.classList.add("uploader");
        document.head.appendChild(script);
    }
}