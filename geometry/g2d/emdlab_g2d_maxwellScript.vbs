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

call defineGlobalVariable(oProject, "x_pts", "[0,36.9154439226,38.4120159735,39.1489850054,51.114594693,52,62.5,60.3703641431,35.7392555727,37,39.6604391634,14,32,36,36.5,32,33.7724994626,34.7135628448,31.6634450549,13.3147912281,15.9371454255] mm")
call makeGBHidden(oProject, "x_pts")
call defineGlobalVariable(oProject, "y_pts", "[0,2.5,2.60135135135,6.34881151725,9.5550096478,0,0,16.1761903189,9.57630466879,0,0,0,0,0,0,9.18794210088,9.18794210088,11.2791202947,10.2880769489,4.32623792125,0] mm")
call makeGBHidden(oProject, "y_pts")
call defineGlobalVariable(oProject, "e_angles", "[0,0,0,-10.5882742142,0,15,0,-11.1257117937,0,0,3.8742882063,9.21149351193,0,0,0,18,0,0,-18,0,1.98,27.2554107456,0] deg")
call makeGBHidden(oProject, "e_angles")
call drawSegment(oEditor, 2, 3, "stator_loop_1_e1")
call drawSegment(oEditor, 3, 4, "stator_loop_1_e2")
call drawSegment(oEditor, 4, 5, "stator_loop_1_e3")
call drawArcCPA(oEditor, 1, 5, 4, "stator_loop_1_e4")
call drawSegment(oEditor, 6, 7, "stator_loop_1_e5")
call drawArcCPA(oEditor, 1, 7, 6, "stator_loop_1_e6")
call drawSegment(oEditor, 8, 9, "stator_loop_1_e7")
call drawArcCPA(oEditor, 1, 9, 8, "stator_loop_1_e8")
call uniteEdges(oEditor, "stator_loop_1_e1,stator_loop_1_e2,stator_loop_1_e3,stator_loop_1_e4,stator_loop_1_e5,stator_loop_1_e6,stator_loop_1_e7,stator_loop_1_e8")
call coverLoop(oEditor, "stator_loop_1_e1")
call rename(oEditor, "stator_loop_1_e1", "stator")
call changeObjectColor(oEditor, "stator", 200,200,200)
call drawSegment(oEditor, 11, 6, "sc_loop_1_e10")
call drawArcCPA(oEditor, 1, 5, 4, "sc_loop_1_e4")
call drawSegment(oEditor, 4, 5, "sc_loop_1_e3")
call drawArcCPA(oEditor, 1, 11, 12, "sc_loop_1_e12")
call uniteEdges(oEditor, "sc_loop_1_e10,sc_loop_1_e4,sc_loop_1_e3,sc_loop_1_e12")
call coverLoop(oEditor, "sc_loop_1_e10")
call rename(oEditor, "sc_loop_1_e10", "sc")
call changeObjectColor(oEditor, "sc", 255,137,39)
call drawSegment(oEditor, 10, 11, "sap_loop_1_e9")
call drawArcCPA(oEditor, 1, 11, 12, "sap_loop_1_e12")
call drawSegment(oEditor, 3, 4, "sap_loop_1_e2")
call drawSegment(oEditor, 2, 3, "sap_loop_1_e1")
call drawArcCPA(oEditor, 1, 10, 11, "sap_loop_1_e11")
call uniteEdges(oEditor, "sap_loop_1_e9,sap_loop_1_e12,sap_loop_1_e2,sap_loop_1_e1,sap_loop_1_e11")
call coverLoop(oEditor, "sap_loop_1_e9")
call rename(oEditor, "sap_loop_1_e9", "sap")
call changeObjectColor(oEditor, "sap", 0,255,255)
call drawSegment(oEditor, 12, 13, "rotor_loop_1_e13")
call drawSegment(oEditor, 13, 16, "rotor_loop_1_e20")
call drawArcCPA(oEditor, 1, 16, 21, "rotor_loop_1_e21")
call drawSegment(oEditor, 19, 20, "rotor_loop_1_e18")
call drawArcCPA(oEditor, 1, 20, 19, "rotor_loop_1_e19")
call uniteEdges(oEditor, "rotor_loop_1_e13,rotor_loop_1_e20,rotor_loop_1_e21,rotor_loop_1_e18,rotor_loop_1_e19")
call coverLoop(oEditor, "rotor_loop_1_e13")
call rename(oEditor, "rotor_loop_1_e13", "rotor")
call changeObjectColor(oEditor, "rotor", 200,200,200)
call drawSegment(oEditor, 13, 14, "magnet_loop_1_e14")
call drawArcCPA(oEditor, 21, 14, 22, "magnet_loop_1_e22")
call drawSegment(oEditor, 17, 16, "magnet_loop_1_e23")
call drawSegment(oEditor, 13, 16, "magnet_loop_1_e20")
call uniteEdges(oEditor, "magnet_loop_1_e14,magnet_loop_1_e22,magnet_loop_1_e23,magnet_loop_1_e20")
call coverLoop(oEditor, "magnet_loop_1_e14")
call rename(oEditor, "magnet_loop_1_e14", "magnet")
call changeObjectColor(oEditor, "magnet", 28,255,28)
call drawSegment(oEditor, 14, 15, "rap_loop_1_e15")
call drawArcCPA(oEditor, 1, 15, 16, "rap_loop_1_e16")
call drawSegment(oEditor, 18, 19, "rap_loop_1_e17")
call drawArcCPA(oEditor, 1, 16, 21, "rap_loop_1_e21")
call drawSegment(oEditor, 17, 16, "rap_loop_1_e23")
call drawArcCPA(oEditor, 21, 14, 22, "rap_loop_1_e22")
call uniteEdges(oEditor, "rap_loop_1_e15,rap_loop_1_e16,rap_loop_1_e17,rap_loop_1_e21,rap_loop_1_e23,rap_loop_1_e22")
call coverLoop(oEditor, "rap_loop_1_e15")
call rename(oEditor, "rap_loop_1_e15", "rap")
call changeObjectColor(oEditor, "rap", 0,255,255)

