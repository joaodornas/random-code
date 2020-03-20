function mySaveAsXML(doc, path)
{
doc.saveAs(path + doc + ".xml","com.adobe.acrobat.xml-1-00");
}
myFunc = app.trustedFunction( function (doc, path)
{
// privileged and/or non-privileged code here
app.beginPriv();
mySaveAsXML(doc, path);
app.endPriv();
// privileged and/or non-privileged code here
});