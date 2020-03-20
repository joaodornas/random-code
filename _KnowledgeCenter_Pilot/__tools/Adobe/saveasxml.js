var inFolder = new Folder("C:/Input")
if(inFolder != null){
var fileList = inFolder.getFiles(/\.(jpg|tif|psd|bmp|gif|png|)$/i);
 }
for(var a = 0 ;a < fileList.length; a++)
{
var docRef = open(fileList[a]);
 //do things here
}

var oRetn = app.browseForDoc({
bSave: true,
cFilenameInit: "myComDoc.pdf",
cFSInit: "CHTTP",
});
if ( typeof oRetn != "undefined" ) this.saveAs({
cFS: oRetn.cFS, cPath: oRetn.cPath, bPromptToOverwrite: false});

var mySaveAs = app.trustedFunction(
   function(oDoc,cPath,cFlName)
   {
      app.beginPriv();
      // Ensure path has trailing "/"
      cPath = cPath.replace(/([^/])$/, "$1/");
      try{
         oDoc.saveAs(cPath + cFlName);
      }catch(e){
         app.alert("Error During Save");
      }
       app.endPriv();
   }
);