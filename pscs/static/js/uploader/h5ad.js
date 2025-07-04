

let h5adUploader = document.getElementById("h5ad");
h5adUploader.addEventListener("change", h5adListener);


function h5adListener(event) {
    const quantEl = document.getElementById("quantities");
    let file = quantEl.files[0];
    if (file) {
        const doTranspose = document.getElementById("quant_transpose").checked;
        let dataNameEl = document.getElementById("file_name");
        if (dataNameEl.value === "") {
            dataNameEl.placeholder = file.name;
        }
    }
}


document.getElementById("form").addEventListener("submit", function (event) {
    document.body.style.cursor = "wait";
    event.target.addEventListener("submit", function () {
        document.body.style.cursor = "default";
    }, {once: true});
});
