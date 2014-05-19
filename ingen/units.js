<?php
if ($GLOBALS['F_UNITS']!="defined") 
{
$GLOBALS['F_UNITS']="defined";

// internally, matrices are "stored" with ps as unit of time and K as unit of energy.

//handy functions to generate menu items
  function option($var,$value,$label,$style="")
  {
     echo "<option value='".$value."' ".($var==$value?"selected='selected'":"")." style='".$style."'>".$label."</option>";
  }

//generate a drop down for units
function unitfield($var,$value,$mode)
{
  echo "<select name='".$value."' onchange='this.form.submit()'> ";
  switch ($mode) {
  case "time":
    option($var,"ps",   "picosec.");
    option($var,"fs",   "femtosec.");
    option($var,"s",    "second");
    option($var,"aut",   "at.time");
    break;
  case "freq":
    option($var,"ps1",   "ps^-1");
    option($var,"thz",   "THz");
    option($var,"thzrad",   "THz rad");
    option($var,"cm1",     "cm^-1");
    option($var,"auw",    "at.freq");
    break;
  case "nrg":
    option($var,"ev",   "eV");
    option($var,"k",    "K");
    option($var,"aue",    "at.energy");
    break;
  }
  echo "</select>";
}

function field($var,$value,$label,$uval,$units)
{
  echo "<div class='label'>".$label."</div>&nbsp;<input class='number' type='text' name='".$value."' value='".$var."' onchange='this.form.submit()'></input>";
  unitfield($uval,"u".$value,$units);
  echo "<br/>\n";
}

function unitname($unit)
{
  switch ($unit) {
  case "ps": return "picoseconds";
  case "fs": return "femtoseconds";
  case "s": return "seconds";
  case "aut": return "atomic time units";
  case "ps1": 
  case "thz": return "THz";
  case "thzrad": return "THz*rad";
  case "cm1": return "cm^-1";
  case "auw": return "atomic frequency units";
  case "k": return "K";
  case "ev": return "eV";
  case "aue": return "atomic energy units"; 
  default: return "Unknown unit $unit";
}
}


function conv_u2i($val,$unit,$mode)
{
  $cc=1.;
  switch ($mode) {
  case "time":
    switch($unit) {
      case "ps": $cc=1.;             break;
      case "fs": $cc=1e-3;             break;
      case "s":  $cc=1.e12;         break; 
      case "aut": $cc=2.4188843e-05;  break;
      default: echo " ERROR! Wrong unit for mode $mode <br/>";
    }
  break;
  case "freq":  
    switch($unit) {
      case "ps1": $cc=1.;             break;
      case "thz": $cc=1.;             break;
      case "thzrad": $cc=1./6.2831853;             break;
      case "cm1":  $cc=0.029979246;         break; 
      case "auw": $cc=6579.6839;  break;
      default: echo " ERROR! Wrong unit for mode $mode <br/>";
    }
  break;
  case "nrg":  
    switch($unit) {
      case "k": $cc=1.;             break;
      case "ev": $cc=11604.506;             break;
      case "kcalmol": $cc=503.55573;             break;
      case "kjmol":  $cc=120.27222;         break; 
      case "aue": $cc=315774.67;  break;
      default: echo " ERROR! Wrong unit for mode $mode <br/>";
    }
  break;
  }
  return $cc*$val;   
}

function conv_i2u($val,$unit,$mode)
{
  $cc=conv_u2i(1.0,$unit,$mode);
  return $val/$cc;
}


}  // ENDS INCLUSION CHECK
?>
