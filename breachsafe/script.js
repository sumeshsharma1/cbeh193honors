  document.addEventListener("mouseover", function(){


/*
chrome.storage.local.get("masterKey", function(data) {
    document.getElementById("demo4").innerHTML = data.masterKey;

});
*/

    chrome.storage.local.get("toStore", function (data) {
    document.getElementById("fields").innerHTML = data.toStore.fields.slice(0,4);
    document.getElementById("time").innerHTML = data.toStore.time;
    document.getElementById("url").innerHTML = data.toStore.url;
 });


});
