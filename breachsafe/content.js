document.addEventListener("submit", function (e) {
  var key = location.host.toString();
  var time = new Date($.now());
  var url = window.location.href;
  var parse = $(e.target).serialize().split("&");
  var info = [];
  for (var i = 0; i < parse.length; i++) {
  		info.push(parse[i].split("=")[0]);
  };

  var toStore = {
  	"fields": info,
  	"time": time.toString(),
  	"url": key
  };

/*
  var masterKey = 42;


var keys = [1,2,3];

chrome.storage.local.get("masterKey", function(data) {
    keys.push([data.masterKey]);
});

keys.push([key + ";" + time.toString()]);

chrome.storage.local.set({masterKey: keys});
*/

chrome.storage.local.set({toStore: toStore});

  chrome.storage.sync.set({[key + "; " + time.toString()]: {toStore}}, function() {
        // Notify that we saved.
        //alert('Submission saved');

		chrome.storage.sync.get([key + "; " + time.toString()], function (data) { console.info(data) });
  });
}, false);


// Remember to change this to the relative path to inject.js
injectScript( chrome.runtime.getURL( "/" ), "inject.js" );

function injectScript ( aBasePath, aScriptURL ) {
  var scriptEl = document.createElement( "script" );
  scriptEl.src = aBasePath + aScriptURL;
  scriptEl.async = false;

  (document.body || document.head || document.documentElement)
  .appendChild( scriptEl );
}
