var pageQuery=function () { 
  var query_string = {};
  var query = window.location.search.substring(1);
  var vars = query.split("&");
  for (var i=0;i<vars.length;i++) {
    var pair = vars[i].split("=");
    	// If first entry with this name
    if (typeof query_string[pair[0]] === "undefined") {
      query_string[pair[0]] = pair[1];
    	// If second entry with this name
    } else if (typeof query_string[pair[0]] === "string") {
      var arr = [ query_string[pair[0]], pair[1] ];
      query_string[pair[0]] = arr;
    	// If third or later entry with this name
    } else {
      query_string[pair[0]].push(pair[1]);
    }
  } 
    return query_string;
} ();

if (!pageQuery.page) pageQuery.page="main"; //default

function httpGet(theUrl, entity="undefined", append=false, silent=true)
{
   /* Reads a file from URL and returns its content as a string*/
   var xmlhttp;
   if (window.XMLHttpRequest)
   {// code for IE7+, Firefox, Chrome, Opera, Safari
      xmlhttp=new XMLHttpRequest();
   }
   else
   {// code for IE6, IE5
      xmlhttp=new ActiveXObject("Microsoft.XMLHTTP");
   }

   if (entity!="undefined") {
      xmlhttp.onreadystatechange=function() {
         if (append) prefix=document.getElementById(entity).innerHTML+" "; else prefix=""; 
         if (xmlhttp.readyState == 4 && (xmlhttp.status == 200 || (xmlhttp.status == 0 && xmlhttp.responseText))) {
            document.getElementById(entity).innerHTML = prefix + xmlhttp.responseText;      
         } else {
            if (!silent) document.getElementById(entity).innerHTML = prefix + "<b> Error. HTTP " + xmlhttp.status + " </b>";           
         }      
      }
   }      
   try {
      xmlhttp.overrideMimeType("text/plain; charset=utf-8");
      xmlhttp.open("GET", theUrl, entity!="undefined" );    
      xmlhttp.send();    
   } catch(err) {  if (!silent) { document.write("Could not load resources from "+theUrl); } throw(err); }

   return xmlhttp.responseText;   
}

function httpInclude(theUrl)
{  /* Reads from URL and writes to document in place */
   document.write(httpGet(theUrl));
}

function composeHeader(pageName)
{
   try { httpGet("heads/base.html", "header", true); }
   catch(err) { document.write("Could not load page header."); }

   try { httpGet("heads/"+pageName+".html", "header", true); }
   catch(err) { httpGet("heads/default.html", "header", true) }   
}

function mkMenu(item,label)
{
   var menu=document.getElementById("menu");
   menu.innerHTML = menu.innerHTML + " <a class='menu " + (item==pageQuery.page?"sel":"")+ "' href='index.html?page="+item+"'>"+label+"</a>";
}

