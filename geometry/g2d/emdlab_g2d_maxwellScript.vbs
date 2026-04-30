sub defineGlobalVariable(oProject, varName, varValue)
oProject.ChangeProperty Array("NAME:AllTabs", Array("NAME:ProjectVariableTab", Array("NAME:PropServers",  _
  "ProjectVariables"), Array("NAME:NewProps", Array("NAME:$"+varName, "PropType:=", "VariableProp", "UserDef:=",  _
  true, "Value:=", varValue))))
end sub

sub makeGBHidden(oProject, varName)
oProject.ChangeProperty Array("NAME:AllTabs", Array("NAME:ProjectVariableTab", Array("NAME:PropServers",  _
  "ProjectVariables"), Array("NAME:ChangedProps", Array("NAME:$"+varName, "Hidden:=", true))))
end sub

sub uniteEdges(oEditor, eNames)

oEditor.Unite Array("NAME:Selections", "Selections:=",  _
  eNames), Array("NAME:UniteParameters", "KeepOriginals:=",  _
  false)

end sub 

sub coverLoop(oEditor, lName)
oEditor.CoverLines Array("NAME:Selections", "Selections:=", lName, "NewPartsModelFlag:=",  _
  "Model")
end sub

sub rename(oEditor, oldName, newName)
oEditor.ChangeProperty Array("NAME:AllTabs", Array("NAME:Geometry3DAttributeTab", Array("NAME:PropServers",  _
  oldName), Array("NAME:ChangedProps", Array("NAME:Name", "Value:=", newName))))
end sub

sub subtract(oEditor, toolParts, blankParts)
oEditor.Subtract Array("NAME:Selections", "Blank Parts:=", blankParts, "Tool Parts:=",  _
  toolParts), Array("NAME:SubtractParameters", "KeepOriginals:=", false)
end sub

sub changeObjectColor(oEditor, oName, R, G, B)
oEditor.ChangeProperty Array("NAME:AllTabs", Array("NAME:Geometry3DAttributeTab", Array("NAME:PropServers",  _
  oName), Array("NAME:ChangedProps", Array("NAME:Color", "R:=", R, "G:=", G, "B:=", B))))
end sub

sub drawSegment(oEditor, index1, index2, name)

oEditor.CreatePolyline Array("NAME:PolylineParameters", "IsPolylineCovered:=", true, "IsPolylineClosed:=",  _
  false, Array("NAME:PolylinePoints", Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index1-1)+"]", "Y:=", "$y_pts["+cstr(index1-1)+"]", "Z:=",  _
  "0mm"), Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index2-1)+"]", "Y:=", "$y_pts["+cstr(index2-1)+"]", "Z:=", "0mm")), Array("NAME:PolylineSegments", Array("NAME:PLSegment", "SegmentType:=",  _
  "Line", "StartIndex:=", 0, "NoOfPoints:=", 2)), Array("NAME:PolylineXSection", "XSectionType:=",  _
  "None", "XSectionOrient:=", "Auto", "XSectionWidth:=", "0mm", "XSectionTopWidth:=",  _
  "0mm", "XSectionHeight:=", "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=",  _
  "Corner")), Array("NAME:Attributes", "Name:=", name, "Flags:=", "", "Color:=",  _
  "(143 175 143)", "Transparency:=", 0, "PartCoordinateSystem:=", "Global", "UDMId:=",  _
  "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SurfaceMaterialValue:=",  _
  "" & Chr(34) & "" & Chr(34) & "", "SolveInside:=", true, "ShellElement:=",  _
  false, "ShellElementThickness:=", "0mm", "IsMaterialEditable:=", true, "UseMaterialAppearance:=",  _
  false, "IsLightweight:=", false)

end sub

sub drawArcCPA(oEditor, index1, index2, index3, name)

oEditor.CreatePolyline Array("NAME:PolylineParameters", "IsPolylineCovered:=", true, "IsPolylineClosed:=",  _
  false, Array("NAME:PolylinePoints", Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index2-1)+"]", "Y:=", "$y_pts["+cstr(index2-1)+"]", "Z:=",  _
  "0mm"), Array("NAME:PLPoint", "X:=", "1mm", "Y:=",  _
  "1mm", "Z:=", "0mm"), Array("NAME:PLPoint", "X:=",  _
  "-0.368087051817776mm", "Y:=", "-0.206384113560369mm", "Z:=", "0mm")), Array("NAME:PolylineSegments", Array("NAME:PLSegment", "SegmentType:=",  _
  "AngularArc", "StartIndex:=", 0, "NoOfPoints:=", 3, "NoOfSegments:=", "0", "ArcAngle:=",  _
  "$e_angles["+cstr(index3-1)+"]", "ArcCenterX:=", "$x_pts["+cstr(index1-1)+"]", "ArcCenterY:=", "$y_pts["+cstr(index1-1)+"]", "ArcCenterZ:=",  _
  "0mm", "ArcPlane:=", "XY")), Array("NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=",  _
  "Auto", "XSectionWidth:=", "0mm", "XSectionTopWidth:=", "0mm", "XSectionHeight:=",  _
  "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner")), Array("NAME:Attributes", "Name:=",  _
  name, "Flags:=", "", "Color:=", "(143 175 143)", "Transparency:=", 0, "PartCoordinateSystem:=",  _
  "Global", "UDMId:=", "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SurfaceMaterialValue:=",  _
  "" & Chr(34) & "" & Chr(34) & "", "SolveInside:=", true, "ShellElement:=",  _
  false, "ShellElementThickness:=", "0mm", "IsMaterialEditable:=", true, "UseMaterialAppearance:=",  _
  false, "IsLightweight:=", false)

end sub

sub drawArc(oEditor, index1, index2, index3, name)

arg1 = "atan2(" + "$y_pts["+cstr(index3-1)+"]" + "-" + "$y_pts["+cstr(index1-1)+"]," + "$x_pts["+cstr(index3-1)+"]" + "-" + "$x_pts["+cstr(index1-1)+"])"
arg2 = "atan2(" + "$y_pts["+cstr(index2-1)+"]" + "-" + "$y_pts["+cstr(index1-1)+"]," + "$x_pts["+cstr(index2-1)+"]" + "-" + "$x_pts["+cstr(index1-1)+"])"

oEditor.CreatePolyline Array("NAME:PolylineParameters", "IsPolylineCovered:=", true, "IsPolylineClosed:=",  _
  false, Array("NAME:PolylinePoints", Array("NAME:PLPoint", "X:=", "$x_pts["+cstr(index2-1)+"]", "Y:=", "$y_pts["+cstr(index2-1)+"]", "Z:=",  _
  "0mm"), Array("NAME:PLPoint", "X:=", "1mm", "Y:=",  _
  "1mm", "Z:=", "0mm"), Array("NAME:PLPoint", "X:=",  _
  "-0.368087051817776mm", "Y:=", "-0.206384113560369mm", "Z:=", "0mm")), Array("NAME:PolylineSegments", Array("NAME:PLSegment", "SegmentType:=",  _
  "AngularArc", "StartIndex:=", 0, "NoOfPoints:=", 3, "NoOfSegments:=", "0", "ArcAngle:=",  _
  arg1+"-"+arg2, "ArcCenterX:=", "$x_pts["+cstr(index1-1)+"]", "ArcCenterY:=", "$y_pts["+cstr(index1-1)+"]", "ArcCenterZ:=",  _
  "0mm", "ArcPlane:=", "XY")), Array("NAME:PolylineXSection", "XSectionType:=", "None", "XSectionOrient:=",  _
  "Auto", "XSectionWidth:=", "0mm", "XSectionTopWidth:=", "0mm", "XSectionHeight:=",  _
  "0mm", "XSectionNumSegments:=", "0", "XSectionBendType:=", "Corner")), Array("NAME:Attributes", "Name:=",  _
  name, "Flags:=", "", "Color:=", "(143 175 143)", "Transparency:=", 0, "PartCoordinateSystem:=",  _
  "Global", "UDMId:=", "", "MaterialValue:=", "" & Chr(34) & "vacuum" & Chr(34) & "", "SurfaceMaterialValue:=",  _
  "" & Chr(34) & "" & Chr(34) & "", "SolveInside:=", true, "ShellElement:=",  _
  false, "ShellElementThickness:=", "0mm", "IsMaterialEditable:=", true, "UseMaterialAppearance:=",  _
  false, "IsLightweight:=", false)

end sub

Dim oAnsoftApp
Dim oDesktop
Dim oProject
Dim oDesign
Dim oEditor
Dim oModule
Set oAnsoftApp = CreateObject("Ansoft.ElectronicsDesktop")
Set oDesktop = oAnsoftApp.GetAppDesktop()
oDesktop.RestoreWindow
Set oProject = oDesktop.NewProject
oProject.InsertDesign "Maxwell 2D", "NewDesign", "Magnetostatic", ""
Set oDesign = oProject.SetActiveDesign("NewDesign")
Set oEditor = oDesign.SetActiveEditor("3D Modeler")

call defineGlobalVariable(oProject, "x_pts", "[43.5,22.4,0,18.786220722,38.1595108415,45.4376438313,20.372042509,16.5039466852,38.4283072328,47.23979003,18.2992028748,14.2994856908,36.9144331692,7.43063708941,2.87839415674,13.6738647888,33.8130895592,9.45998193221,4.24139849268,14.4249783362,30.7591449816,2.66360678815e-15,1.37160441505e-15,-12.1999143843,-24.7810760932,2.38128600905,1.06765350736,-11.9908191468,-27.9197994839,4.96510200521,1.92332372412,-11.5794951953,-29.8927185749,-46.9151961783,-18.173465467,-12.312003157,-30.4454425893,-44.5057158334,-19.954211055,-14.4249783362,-30.7591449816,-43.5,-22.4,-18.786220722,-38.1595108415,-45.4376438313,-20.372042509,-16.5039466852,-38.4283072328,-47.23979003,-18.2992028748,-14.2994856908,-36.9144331692,-7.43063708941,-2.87839415674,-13.6738647888,-33.8130895592,-9.45998193221,-4.24139849268,-14.4249783362,-30.7591449816,-7.99082036444e-15,-4.11481324514e-15,12.1999143843,24.7810760932,-2.38128600905,-1.06765350736,11.9908191468,27.9197994839,-4.96510200521,-1.92332372412,11.5794951953,29.8927185749,46.9151961783,18.173465467,12.312003157,30.4454425893,44.5057158334,19.954211055,14.4249783362,30.7591449816] mm")
call makeGBHidden(oProject, "x_pts")
call defineGlobalVariable(oProject, "y_pts", "[0,0,0,12.1999143843,24.7810760932,-2.38128600905,-1.06765350736,11.9908191468,27.9197994839,-4.96510200521,-1.92332372412,11.5794951953,29.8927185749,46.9151961783,18.173465467,12.312003157,30.4454425893,44.5057158334,19.954211055,14.4249783362,30.7591449816,43.5,22.4,18.786220722,38.1595108415,45.4376438313,20.372042509,16.5039466852,38.4283072328,47.23979003,18.2992028748,14.2994856908,36.9144331692,7.43063708941,2.87839415674,13.6738647888,33.8130895592,9.45998193221,4.24139849268,14.4249783362,30.7591449816,5.32721357629e-15,2.74320883009e-15,-12.1999143843,-24.7810760932,2.38128600905,1.06765350736,-11.9908191468,-27.9197994839,4.96510200521,1.92332372412,-11.5794951953,-29.8927185749,-46.9151961783,-18.173465467,-12.312003157,-30.4454425893,-44.5057158334,-19.954211055,-14.4249783362,-30.7591449816,-43.5,-22.4,-18.786220722,-38.1595108415,-45.4376438313,-20.372042509,-16.5039466852,-38.4283072328,-47.23979003,-18.2992028748,-14.2994856908,-36.9144331692,-7.43063708941,-2.87839415674,-13.6738647888,-33.8130895592,-9.45998193221,-4.24139849268,-14.4249783362,-30.7591449816] mm")
call makeGBHidden(oProject, "y_pts")
call defineGlobalVariable(oProject, "e_angles", "[0,33,0,-36,0,39,0,-42,0,45,0,42,0,-39,0,36,0,-33,0,0,33,0,-36,0,39,0,-42,0,45,0,42,0,-39,0,36,0,-33,0,0,33,0,-36,0,39,0,-42,0,45,0,42,0,-39,0,36,0,-33,0,0,33,0,-36,0,39,0,-42,0,45,0,42,0,-39,0,36,0,-33,0] deg")
call makeGBHidden(oProject, "e_angles")

